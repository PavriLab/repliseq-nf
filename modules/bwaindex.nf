process BWA_INDEX {

    tag "$bwa_base"

    input:
    file(fasta)

    output:
    path("BWAIndex")

    shell:
    bwa_base = fasta.toString() - ~/(\.fa)?(\.fasta)?(\.fas)?$/

    """
    mkdir bowtie2Index

    bwa index -a bwtsw ${fasta}
    mkdir BWAIndex && mv ${fasta}* BWAIndex
    """
}
