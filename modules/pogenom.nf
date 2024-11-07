process POGENOM {

    conda "/home/daniel/.conda/envs/pogenom"

    publishDir "${params.outdir}/pogenom", mode: 'copy'

    input:
    tuple val(sample_id), path(mag), path(bam_file), path(vcf_file), path(coverage_file)

    output:
    tuple val(sample_id), path("${sample_id}_pogenom_results"), emit: results

    script:
    """
    pogenom --genome_file ${mag} --vcf_file ${vcf_file} --coverage_file ${coverage_file} --output_dir ${sample_id}_${mag.baseName}_pogenom_results
    """
}
