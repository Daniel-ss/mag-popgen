process BOWTIE2_BUILD {
    conda "/opt/miniforge3/envs/mapeo_cacao"
    
    publishDir "${params.outdir}/bowtie2/build", mode: 'copy'

    input:
    tuple val(sample_id), path(mag)
     
    output:
    tuple val(sample_id), path("${sample_id}_${mag.baseName}*"), emit: index

    script:
    """
    ## Build index
    bowtie2-build \
     --threads 12 \
     ${mag} \
     ${sample_id}_${mag.baseName}

    """
}
