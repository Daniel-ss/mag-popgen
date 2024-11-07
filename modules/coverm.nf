process COVERM {
    conda "/opt/miniforge3/envs/genomics"

    publishDir "${params.outdir}/coverm", mode: 'copy'

    input:
    tuple val(reference_id), path(mag), path(bam_files), path(contig_names)

    output:
    tuple val(reference_id), path(mag), path("${reference_id}_coverage.txt"), emit: coverage

    script:
    """
    coverm genome \
     --genome-definition ${contig_names} \
     --bam-files ${bam_files.join(' ')} \
     --methods relative_abundance \
     --output-file ${reference_id}_coverage.txt
    """
}
