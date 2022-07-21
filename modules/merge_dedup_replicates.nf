process MERGE_DEDUP_REPLICATES {

    tag "$meta.id"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filemeta.id ->
                      if (filemeta.id.endsWith(".flagstat")) "samtools_stats/$filemeta.id"
                      else if (filemeta.id.endsWith(".idxstats")) "samtools_stats/$filemeta.id"
                      else if (filemeta.id.endsWith(".stats")) "samtools_stats/$filemeta.id"
                      else if (filemeta.id.endsWith(".metrics.txt")) "picard_metrics/$filemeta.id"
                      else if (filemeta.id.endsWith(".bam")) "bams/$filemeta.id"
                      else if (filemeta.id.endsWith(".bai")) "bams/$filemeta.id"
                      else filemeta.id
                }

    input:
    tuple val(meta), file(alignments)

    output:
    tuple val(meta), file("*.markdup.{bam,bam.bai}"),   emit: alignments
    tuple val(meta), file("*.flagstat"),                emit: flagstats
    path("*.{idxstats,stats}"),                         emit: bamstats
    path("*.txt"),                                      emit: dedupstats

    script:
    bam_files = alignments.findAll { it.toString().endsWith('.bam') }.sort()

    if (bam_files.size() > 1) {
        """
        samtools merge -@ $task.cpus ${meta.id}.merged.bam ${bam_files.join(' ')}

        samtools index ${meta.id}.merged.bam

        picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
            INPUT=${meta.id}.merged.bam \\
            OUTPUT=${meta.id}.markdup.bam \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=true \\
            METRICS_FILE=${meta.id}.MarkDuplicates.metrics.txt \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp

        samtools index ${meta.id}.markdup.bam
        samtools flagstat ${meta.id}.markdup.bam > ${meta.id}.markdup.bam.flagstat
        samtools idxstats ${meta.id}.markdup.bam > ${meta.id}.markdup.bam.idxstats
        samtools stats ${meta.id}.markdup.bam > ${meta.id}.markdup.bam.stats
        """
    } else {
      """

      picard -Xmx${task.memory.toGiga()}g MarkDuplicates \\
          INPUT=${alignments[0]} \\
          OUTPUT=${meta.id}.markdup.bam \\
          ASSUME_SORTED=true \\
          REMOVE_DUPLICATES=true \\
          METRICS_FILE=${meta.id}.MarkDuplicates.metrics.txt \\
          VALIDATION_STRINGENCY=LENIENT \\
          TMP_DIR=tmp

      samtools index ${meta.id}.markdup.bam
      samtools flagstat ${meta.id}.markdup.bam > ${meta.id}.markdup.bam.flagstat
      samtools idxstats ${meta.id}.markdup.bam > ${meta.id}.markdup.bam.idxstats
      samtools stats ${meta.id}.markdup.bam > ${meta.id}.markdup.bam.stats
      """
    }
}
