process CHECKM2 {

    conda "/opt/miniforge3/envs/checkm2"

    publishDir "${params.outdir}/checkm2", mode: 'copy'

    input:
    tuple val(sample_id), path(mag)

    output:
    tuple val(sample_id), path("${sample_id}_${mag.baseName}_checkm2_results"), emit: results

    script:
    """

    checkm2 predict \
     --threads ${task.cpus} \
     --input ${mag} \
     --output-directory ${sample_id}_${mag.baseName}_checkm2_results \
     --database_path ${params.checkm2_db}
    """
}
