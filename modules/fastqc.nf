process FASTQC {

   tag "$meta.id"

   input:
   tuple val(meta), file(reads)

   output:
   tuple val(meta), file("*.{zip,html}"), emit: multiqc

   script:
   if (params.singleEnd) {
       """
       [ ! -f  ${meta.id}.fastq.gz ] && ln -s $reads ${meta.id}.fastq.gz
       fastqc -q -t $task.cpus ${meta.id}.fastq.gz
       """

   } else {
       """
       [ ! -f  ${meta.id}_1.fastq.gz ] && ln -s ${reads[0]} ${meta.id}_1.fastq.gz
       [ ! -f  ${meta.id}_2.fastq.gz ] && ln -s ${reads[1]} ${meta.id}_2.fastq.gz
       fastqc -q -t $task.cpus ${meta.id}_1.fastq.gz
       fastqc -q -t $task.cpus ${meta.id}_2.fastq.gz
       """
   }
}
