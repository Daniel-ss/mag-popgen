process PROKKA {
        
    conda "/opt/miniforge3/envs/genomics"

    publishDir "${params.outdir}/prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(mag)

    output:
    tuple val(sample_id), path("${sample_id}_${mag.baseName}_annotation"), emit: annotation

    script:
    """
    prokka --outdir ${sample_id}_${mag.baseName}_annotation --prefix ${sample_id}_${mag.baseName} ${mag}
    """
}
