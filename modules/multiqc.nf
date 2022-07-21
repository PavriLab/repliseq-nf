process MULTIQC {

    tag "multiqc"

    publishDir "${params.outputDir}/multiqc/", mode: 'copy'

    input:
    file(multiqcConfig)

    file ('fastqc/*') from fastqcMultiqcChannel.collect().ifEmpty([])
    file ('trimgalore/*') from trimgaloreMultiqcChannel.collect().ifEmpty([])
    file ('trimgalore/fastqc/*') from trimgaloreFastqcMultiqcChannel.collect().ifEmpty([])

    file ('alignment/replicates/*') from bwaMultiqcChannel.collect()
    file ('alignment/*') from mergeFlagstatMultiqcChannel.collect()
    file ('alignment/*') from mergeidxStatsMultiqcChannel.collect()
    file ('alignment/picard_metrics/*') from mergeMarkDuplicatesMultiqcChannel.collect()

    output:
    path("*multiqc_report.html" into multiQCChannel
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqcConfig \\
        -m custom_content -m fastqc -m cutadapt -m samtools -m picard
    """
}
