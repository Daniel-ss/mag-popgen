process FREEBAYES {
    conda "/opt/miniforge3/envs/genomics"

    publishDir "${params.outdir}/freebayes", mode: 'copy'

    input:
    tuple val(reference_mag), path(mag), path(bams)

    output:
    tuple val(reference_mag), path("${reference_mag}_unfiltered.vcf"), emit: vcf

    script:
    """
    # Index the reference genome
    samtools faidx ${mag}

    # Run freebayes with all BAM files
    freebayes \
     -f ${mag} \
     -C 4 \
     -p 1 \
     -z 0.05 \
     -F 0.01 \
     -q 15 \
     --report-monomorphic \
     ${bams.collect{"-b $it"}.join(' ')} > ${reference_mag}_unfiltered.vcf
    """
}
