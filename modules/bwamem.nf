process BWAMem {

    tag "$meta.id"

    publishDir "${params.outdir}", mode: 'copy',
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
    tuple val(meta), file(reads)
    path(index)
    path(bamtoolsFilterConfig)

    output:
    tuple val(meta), file("*.filtered.{bam,bam.bai}"),  emit: alignments
    path("*.{flagstat,idxstats,stats}"),                emit: multiqc

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
            -out ${meta.id}.filtered.bam \\
            -script $bamtoolsFilterConfig

    samtools index ${meta.id}.filtered.bam
    samtools flagstat ${meta.id}.filtered.bam > ${meta.id}.filtered.bam.flagstat
    samtools idxstats ${meta.id}.filtered.bam > ${meta.id}.filtered.bam.idxstats
    samtools stats ${meta.id}.filtered.bam > ${meta.id}.filtered.bam.stats
    """
}
