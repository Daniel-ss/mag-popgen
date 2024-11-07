process FREEBAYES {

    conda "/opt/miniforge3/envs/genomics"

    publishDir "${params.outdir}/freebayes", mode: 'copy'

    input:
    tuple val(reference_id), path(mag), path(bam_files)

    output:
    tuple val(reference_id), path("${reference_id}_variants.vcf"), emit: vcf

    script:
    """
    # Index the reference genome
    samtools faidx ${mag}

    # Run freebayes with all BAM files
    freebayes \
     -f ${mag} ${bam_files.join(' ')} > \
     ${reference_id}_variants.vcf
    """
}
