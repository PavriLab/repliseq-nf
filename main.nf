#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2020 Tobias Neumann
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

def helpMessage() {
    log.info"""
    ================================================================
     hicer-nf
    ================================================================
     DESCRIPTION

     Basic processing of RepliSeq data.

     Usage:
     nextflow run t-neumann/repliseq-nf

     Options:
        --samples        Tab-delimited text file specifying the samples
                         to be processed. (default: 'samples.txt')
                         The following columns are required:
                            - name: name of sample
                            - read1: Read file with first read mates (R1) in fastq(.gz) format
                            - read2: Read file with second read mates (R2) in fastq(.gz) format

        --singleEnd      Specifies that the input is single-end reads

        --outputDir      Directory name to save results to. (Defaults to
                         'results')

        References:
        --genome         Name of reference (hg38, mm10)
        --fasta          Alternatively, path to genome fasta file which will be digested

     Profiles:
        standard         local execution
        singularity      local execution with singularity
        cbe              CBE cluster execution with singularity

     Docker:
     zuberlab/repliseq-nf:latest

     Authors:
     Tobias Neumann (tobias.neumann@imp.ac.at)
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

params.hicdigest = params.genome ? params.genomes[ params.genome ].hicdigest ?: false : false
params.hicRestriction = params.genome ? params.genomes[ params.genome ].hicRestriction ?: false : false
params.bowtie2 = params.genome ? params.genomes[ params.genome ].bowtie2 ?: false : false

log.info ""
log.info " parameters "
log.info " ======================"
log.info " Design                   : ${params.design}"
log.info " Single end               : ${params.singleEnd}"
log.info " Genome                   : ${params.genome}"
log.info " Fasta                    : ${params.fasta}"
log.info " BWA index                : ${params.bwa}"
log.info " Output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

if (params.design)     { ch_input = file(params.design, checkIfExists: true) } else { exit 1, "Samples design file not specified!" }

if (params.bwa) {

  lastPath = params.bwa.lastIndexOf(File.separator)
  bwa_dir =  params.bwa.substring(0,lastPath+1)
  bwa_base = params.bwa.substring(lastPath+1)

  bwaIndex = Channel.fromPath( bwa_dir , checkIfExists: true)
      .ifEmpty { exit 1, "Genome index: Provided index not found: ${params.bwa}" }

} else if (params.fasta) {
  lastPath = params.fasta.lastIndexOf(File.separator)
  bwa_base = params.fasta.substring(lastPath+1)

  fastaForBwa = Channel.fromPath(params.fasta, checkIfExists: true)
      .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
} else {
    exit 1, "No reference genome files specified!"
}

if (!params.bwa && params.fasta) {
      process BWAIndex {
          tag "$bwa_base"

          input:
          file fasta from fastaForBwa

          output:
          file "BWAIndex" into bwaIndex

          shell:
          bwa_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/
          """
          mkdir bowtie2Index

          bwa index -a bwtsw ${fasta}
          mkdir BWAIndex && mv ${fasta}* BWAIndex

	        """

      }
}

Channel
    .fromPath( params.samples )
    .splitCsv(sep: '\t', header: true)
    .set { samples }

/*
 * PREPROCESSING - REFORMAT DESIGN FILE AND CHECK VALIDITY
 */
process CheckDesign {
    tag "$design"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file design from ch_input

    output:
    file "*.csv" into ch_design_reads_csv

    script:
    """
    check_design.py $design design_reads.csv
    """
}

/*
 * Create channels for input fastq files
 */
if (params.single_end) {
    ch_design_reads_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ] ] }
        .into { ch_raw_reads_fastqc;
                ch_raw_reads_trimgalore;
                design_replicates_exist;
                design_multiple_samples }
} else {
    ch_design_reads_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
        .into { ch_raw_reads_fastqc;
                ch_raw_reads_trimgalore;
                design_replicates_exist;
                design_multiple_samples }
}

// Boolean value for replicates existing in design
replicatesExist = design_replicates_exist
                      .map { it -> it[0].split('_')[-2].replaceAll('R','').toInteger() }
                      .flatten()
                      .max()
                      .val > 1

// Boolean value for multiple groups existing in design
multipleGroups = design_multiple_samples
                     .map { it -> it[0].split('_')[0..-3].join('_') }
                     .flatten()
                     .unique()
                     .count()
                     .val > 1

/*
* STEP 1 - FastQC
*/
process FastQC {
   tag "$name"

   input:
   set val(name), file(reads) from ch_raw_reads_fastqc

   output:
   file "*.{zip,html}" into ch_fastqc_reports_mqc

   script:
   if (params.single_end) {
       """
       [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
       fastqc -q -t $task.cpus ${name}.fastq.gz
       """
   } else {
       """
       [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
       [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
       fastqc -q -t $task.cpus ${name}_1.fastq.gz
       fastqc -q -t $task.cpus ${name}_2.fastq.gz
       """
   }
}

/*
 * STEP 2 - Trim Galore!
 */
process TrimGalore {
    tag "$name"

    input:
    set val(name), file(reads) from ch_raw_reads_trimgalore

    output:
    set val(name), file("*.fq.gz") into ch_trimmed_reads
    file "*.txt" into ch_trimgalore_results_mqc
    file "*.{zip,html}" into ch_trimgalore_fastqc_reports_mqc

    script:
    if (params.single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        trim_galore --fastqc --gzip ${name}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        trim_galore --paired --fastqc --gzip ${name}_1.fastq.gz ${name}_2.fastq.gz
        """
    }
}

/*
 * STEP 3 - Align read with bwa
 */
process BWAMem {
    tag "$name"

    input:
    set val(name), file(reads) from ch_trimmed_reads
    file index from ch_bwa_index.collect()

    output:
    set val(name), file("*.bam") into ch_bwa_bam

    script:
    """
    bwa mem \\
        -t $task.cpus \\
        -M \\
        ${index}/${bwa_base} \\
        $reads \\
        | samtools view -@ $task.cpus -b -h -F 0x0100 -O BAM -o ${name}.bam -
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
