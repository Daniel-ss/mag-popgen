process COVERM {

    publishDir "${params.outdir}/coverm", mode: 'copy'

    input:
    tuple val(sample_id), path(genome_definition), val(bam_files)

    output:
    tuple val(sample_id), path("${sample_id}_coverage.txt"), emit: coverage

    script:
    """
    coverm genome \
     --genome-definition ${genome_definition} \
     --bam-files ${bam_files.join(' ')} \
     --methods relative_abundance \
     --threads 12 \
     --output-file ${sample_id}_coverage.txt
    """
}
