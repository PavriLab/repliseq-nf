process COMPUTE_EL_RATIO {

    tag "$meta.id"

    input:
    tuple val(meta), file(alignments)

    output:
    tuple val(meta) file("*.bg"),   emit: bedgraph

    script:

    earlyBam = alignments[0].containsKey("E") ? alignments[0]["E"] : alignments[1]["E"]
    lateBam  = alignments[0].containsKey("L") ? alignments[0]["L"] : alignments[1]["L"]

    """
    bamCompare -b1 ${earlyBam[0]} \
               -b2 ${lateBam[0]} \
               -o ${meta.id}.bg -of bedgraph \
               -bs ${params.windowSize} \
               --scaleFactorsMethod readCount \
               --operation log2 -p $task.cpus
    """
}
