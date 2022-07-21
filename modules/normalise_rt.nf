process NORMALISE_RT {

    tag "$meta.id"

    publishDir "${params.outputDir}", mode: 'copy', overwrite: 'true',
      saveAs: { filename ->
                    if (filename.endsWith(".bedGraph")) "bedGraph/$filename"
                    else filename
              }

    input:
    tuple val(meta), file(bedgraph)

    output:
    tuple val(meta), file("*.bedGraph"),    emit: bedgraph

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
