process CHROMSIZES_FROM_FASTA {

    tag "$fasta"

    input:
    file(fasta)

    output:

    path("*.sizes")

    script:
    """
    samtools faidx $fasta
    cut -f 1,2 ${fasta}.fai > ${fasta}.sizes
    """
}
