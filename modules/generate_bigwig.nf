process BIGWIG {

    tag "$meta.id"

    publishDir "${params.outputDir}", mode: 'copy', overwrite: 'true',
      saveAs: { filename ->
                    if (filename.endsWith(".qnorm.bw")) "bigwigs/quantile_normalized/$filename"
                    else if (filename.endsWith(".qnorm.loess.bw")) "bigwigs/quantile_normalized_loess/$filename"
                    else if (filename.endsWith(".loess.bw")) "bigwigs/unnormalized_loess/$filename"
                    else "bigwigs/unnormalized/$filename"
              }

    input:
    tuple val(meta), file(bedgraph)
    file(chrsizes)

    output:
    path("*.bw"),   emit: bigwig

    script:
    """
    bedGraphToBigWig $bedgraph $chrsizes ${meta.id}.bw
    """
}
