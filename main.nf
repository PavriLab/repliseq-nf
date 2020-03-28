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
        --design         Tab-delimited text file specifying the samples
                         to be processed. (default: 'samples.txt')
                         The following columns are required:
                            - name: name of sample
                            - read1: Read file with first read mates (R1) in fastq(.gz) format
                            - read2: Read file with second read mates (R2) in fastq(.gz) format

        --singleEnd      Specifies that the input is single-end reads

        --windowSize     Specifies the non-overlapping window size for binning

        --loessSpan      Specifies the span for Loess smoothing

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



log.info ""
log.info " parameters "
log.info " ======================"
log.info " Design                   : ${params.design}"
log.info " Window size              : ${params.windowSize}"
log.info " Loess span               : ${params.loessSpan}"
log.info " Single end               : ${params.singleEnd}"
log.info " Genome                   : ${params.genome}"
log.info " Fasta                    : ${params.fasta}"
log.info " BWA index                : ${params.bwa}"
log.info " Output directory         : ${params.outputDir}"
log.info " ======================"
log.info ""

params.bwa = params.genome ? params.genomes[ params.genome ].bwa ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false

if (params.design)     { designChannel = file(params.design, checkIfExists: true) } else { exit 1, "Samples design file not specified!" }

if (params.singleEnd) {
    bamtoolsFilterConfigChannel = file(params.bamtoolsFilterSEConf, checkIfExists: true)
} else {
    bamtoolsFilterConfigChannel = file(params.bamtoolsFilterPEConf, checkIfExists: true)
}

if (!params.fasta) {
   exit 1, "No genome fasta specified!"
} else {
   fastaForGenomeSizes = Channel.fromPath( params.fasta , checkIfExists: true)
      .ifEmpty { exit 1, "Genome fasta file not found: ${params.fasta}" }
}

multiqcConfigChannel = file(params.multiqcConfig, checkIfExists: true)

customRunName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    customRunName = workflow.runName
}

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

process MakeGenomeFilter {
    tag "$fasta"

    input:
    file fasta from fastaForGenomeSizes

    output:

    file "*.sizes" into chromSizesChannel

    script:
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    """
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

/*
 * PREPROCESSING - REFORMAT DESIGN FILE AND CHECK VALIDITY
 */
process CheckDesign {
    tag "$design"
    publishDir "${params.outputDir}/pipeline_info", mode: 'copy'

    input:
    file design from designChannel

    output:
    file "*.csv" into designCheckedChannel

    script:
    """
    check_design.py $design design_reads.csv
    """
}

/*
 * Create channels for input fastq files
 */
if (params.singleEnd) {
    designCheckedChannel
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ] ] }
        .into { rawReadsFastqcChannel;
                rawReadsTrimgaloreChannel;
                designMultipleGroups }
} else {
    designCheckedChannel
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ] }
        .into { rawReadsFastqcChannel;
                rawReadsTrimgaloreChannel;
                designMultipleGroups }
}

// Boolean value for multiple groups existing in design
multipleGroups = designMultipleGroups
                     .map { it -> it[0].split('_')[0] }
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
   set val(name), file(reads) from rawReadsFastqcChannel

   output:
   file "*.{zip,html}" into fastqcMultiqcChannel

   script:
   if (params.singleEnd) {
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
    set val(name), file(reads) from rawReadsTrimgaloreChannel

    output:
    set val(name), file("*.fq.gz") into trimmedReadChannel
    file "*.txt" into trimgaloreMultiqcChannel
    file "*.{zip,html}" into trimgaloreFastqcMultiqcChannel

    script:
    if (params.singleEnd) {
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
 * STEP 3 - Align reads with bwa and filter
 */
process BWAMem {
    tag "$name"

    input:
    set val(name), file(reads) from trimmedReadChannel
    file index from bwaIndex.collect()
    file bamtoolsFilterConfig from bamtoolsFilterConfigChannel

    output:
    set val(name), file("*.filtered.{bam,bam.bai}") into bwaChannel
    file "*.{flagstat,idxstats,stats}" into bwaMultiqcChannel

    script:
    filter_params = params.singleEnd ? "-F 0x004" : "-F 0x004 -F 0x0008 -f 0x001"
    """
    bwa mem \\
        -t $task.cpus \\
        -M \\
        ${index}/${bwa_base} \\
        $reads \\
        | samtools view -@ $task.cpus -b -h -F 0x0100 -O BAM -o unsorted.bam -

    samtools sort -@ $task.cpus -o sorted.bam unsorted.bam

    samtools view \\
        $filter_params \\
        -F 0x0400 -q 1 \\
        -b sorted.bam \\
        | bamtools filter \\
            -out ${name}.filtered.bam \\
            -script $bamtoolsFilterConfig

    samtools index ${name}.filtered.bam
    samtools flagstat ${name}.filtered.bam > ${name}.filtered.bam.flagstat
    samtools idxstats ${name}.filtered.bam > ${name}.filtered.bam.idxstats
    samtools stats ${name}.filtered.bam > ${name}.filtered.bam.stats

    """
}

/*
 * STEP 7 Merge library BAM files across all replicates
 */
bwaChannel
    .map { it -> [ it[0].split('_')[0] + "_" + it[0].split('_')[1], it[1] ] }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { mergeChannel }

process MergedRepBAM {
    tag "$name"
    publishDir "${params.outputDir}", mode: 'copy',
        saveAs: { filename ->
                      if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                      else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                      else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                      else if (filename.endsWith(".metrics.txt")) "picard_metrics/$filename"
                      else if (filename.endsWith(".bam")) "bams/$filename"
                      else if (filename.endsWith(".bai")) "bams/$filename"
                      else filename
                }

    input:
    set val(name), file(bams) from mergeChannel

    output:
    set val(name), file("*.markdup.{bam,bam.bai}") into bedGraphChannel
    set val(name), file("*.flagstat") into mergeFlagstatMultiqcChannel
    file "*.{idxstats,stats}" into mergeidxStatsMultiqcChannel
    file "*.txt" into mergeMarkDuplicatesMultiqcChannel

    script:
    bam_files = bams.findAll { it.toString().endsWith('.bam') }.sort()

    if (bam_files.size() > 1) {
        """
        samtools merge -@ $task.cpus ${name}.merged.bam ${bam_files.join(' ')}

        samtools index ${name}.merged.bam

        picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
            INPUT=${name}.merged.bam \\
            OUTPUT=${name}.markdup.bam \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=true \\
            METRICS_FILE=${name}.MarkDuplicates.metrics.txt \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp

        samtools index ${name}.markdup.bam
        samtools flagstat ${name}.markdup.bam > ${name}.markdup.bam.flagstat
        samtools idxstats ${name}.markdup.bam > ${name}.markdup.bam.idxstats
        samtools stats ${name}.markdup.bam > ${name}.markdup.bam.stats
        """
    } else {
      """

      picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
          INPUT=${bams[0]} \\
          OUTPUT=${name}.markdup.bam \\
          ASSUME_SORTED=true \\
          REMOVE_DUPLICATES=true \\
          METRICS_FILE=${name}.MarkDuplicates.metrics.txt \\
          VALIDATION_STRINGENCY=LENIENT \\
          TMP_DIR=tmp

      samtools index ${name}.markdup.bam
      samtools flagstat ${name}.markdup.bam > ${name}.markdup.bam.flagstat
      samtools idxstats ${name}.markdup.bam > ${name}.markdup.bam.idxstats
      samtools stats ${name}.markdup.bam > ${name}.markdup.bam.stats
      """
    }
}

bedGraphChannel
    .map { it ->
        def condition = it[0].split('_')[0]
        def phase = it[0].split('_')[1]
        def dictionary = [ (phase) : it[1].flatten()]
        return tuple(condition, dictionary)
     }
    .groupTuple(by: [0])
    .map { it ->  [ it[0], it[1].flatten() ] }
    .set { conditionChannel }

process ELRatio {

    tag "$name"

    input:
    set val(name), val(bam) from conditionChannel

    output:
    file("*.bg") into ELRatioChannel

    script:

    earlyBam = bam[0].containsKey("E") ? bam[0]["E"] : bam[1]["E"]
    lateBam  = bam[0].containsKey("L") ? bam[0]["L"] : bam[1]["L"]

    """
    bamCompare -b1 ${earlyBam[0]} \
               -b2 ${lateBam[0]} \
               -o ${name}.bg -of bedgraph \
               -bs ${params.windowSize} \
               --scaleFactorsMethod readCount \
               --operation log2 -p $task.cpus
    """
}

process RTNormalization {

    tag "$name"

    input:
    file(bedgraph) from ELRatioChannel.collect()

    output:
    file("*.bedGraph") into RTNormalizationChannel

    shell:

    if (multipleGroups) {

      '''
      echo -e "chr\tstart\tstop\t"`ls *.bg` | sed "s/\\ /\\t/g" > merge_RT.txt
      bedtools unionbedg -filler "NA" -i *.bg >> merge_RT.txt

      rtnormalize.r -r merge_RT.txt -s !{params.loessSpan}

      for file in `ls *bg`;
      do
        name=${file%.bg}
        sort -k1,1 -k2,2n $file > ${name}.bedGraph
      done

      '''

    } else {

      '''
      echo -e "chr\tstart\tstop\t"`ls *.bg` | sed "s/\\ /\\t/g" > merge_RT.txt
      cat *.bg >> merge_RT.txt

      rtnormalize.r -r merge_RT.txt -s !{params.loessSpan}

      for file in `ls *bg`;
      do
        name=${file%.bg}
        sort -k1,1 -k2,2n $file > ${name}.bedGraph
      done

      '''
    }
}

RTNormalizationChannel
  .flatten()
  .map { it ->  [ it.baseName, it ] }
  .set { bigwigInputChannel }

process bigwig {

  publishDir "${params.outputDir}", mode: 'copy', overwrite: 'true',
      saveAs: { filename ->
                    if (filename.endsWith(".qnorm.bw")) "bigwigs/quantile_normalized/$filename"
                    else if (filename.endsWith(".qnorm.loess.bw")) "bigwigs/quantile_normalized_loess/$filename"
                    else if (filename.endsWith(".loess.bw")) "bigwigs/unnormalized_loess/$filename"
                    else "bigwigs/unnormalized/$filename"
              }

    tag "$name"

    input:
    file(chrsizes) from chromSizesChannel.collect()
    set val(name), file(bedgraph) from bigwigInputChannel

    output:
    file("*.bw") into bigwigOutput

    script:

    """
    bedGraphToBigWig $bedgraph $chrsizes ${name}.bw
    """
}

process MultiQC {
    publishDir "${params.outputDir}/multiqc/", mode: 'copy'

    input:
    file multiqcConfig from multiqcConfigChannel

    file ('fastqc/*') from fastqcMultiqcChannel.collect().ifEmpty([])
    file ('trimgalore/*') from trimgaloreMultiqcChannel.collect().ifEmpty([])
    file ('trimgalore/fastqc/*') from trimgaloreFastqcMultiqcChannel.collect().ifEmpty([])

    file ('alignment/replicates/*') from bwaMultiqcChannel.collect()
    file ('alignment/*') from mergeFlagstatMultiqcChannel.collect()
    file ('alignment/*') from mergeidxStatsMultiqcChannel.collect()
    file ('alignment/picard_metrics/*') from mergeMarkDuplicatesMultiqcChannel.collect()

    output:
    file "*multiqc_report.html" into multiQCChannel
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename \\
        -m custom_content -m fastqc -m cutadapt -m samtools -m picard
    """
}


workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
