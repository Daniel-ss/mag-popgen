process BOWTIE2_BUILD {
    conda "/opt/miniforge3/envs/mapeo_cacao"
    
    publishDir "${params.outdir}/bowtie2/build", mode: 'copy'

    input:
    tuple val(reference_id), val(mag_id), path(mag)
     
    output:
    tuple val(reference_id), val("${reference_id}_${mag.baseName}"), path("${reference_id}_${mag.baseName}*.bt2"), emit: index

    script:
    """
    ## Build index
    bowtie2-build \
     --threads 12 \
     ${mag} \
     ${reference_id}_${mag.baseName}

    """
}
