process TRIM_GALORE {

    tag "$meta.id"

    input:
    tuple val(meta), file(reads)

    output:
    tuple val(meta), file("*.fq.gz"),   emit: reads
    path("*.txt"),                      emit: multiqc
    path("*.{zip,html}"),               emit: fastqc

    script:
    if (params.singleEnd) {
        """
        [ ! -f  ${meta.id}.fastq.gz ] && ln -s $reads ${meta.id}.fastq.gz
        trim_galore --fastqc --gzip ${meta.id}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${meta.id}_1.fastq.gz ] && ln -s ${reads[0]} ${meta.id}_1.fastq.gz
        [ ! -f  ${meta.id}_2.fastq.gz ] && ln -s ${reads[1]} ${meta.id}_2.fastq.gz
        trim_galore --paired --fastqc --gzip ${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz
        """
    }
}
