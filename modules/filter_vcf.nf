process FILTER_VCF {
    
    publishDir "${params.outdir}/filtered_vcf", mode: 'copy'

    input:
    tuple val(reference_mag), path(vcf)

    output:
    tuple val(reference_mag), path("${reference_mag}_filtered.vcf"), emit: filtered_vcf

    script:
    """
    bcftools filter -i 'QUAL >= 20' ${vcf} > ${reference_mag}_filtered.vcf
    #| bgzip > ${reference_mag}_filtered.vcf.gz
    #tabix -p vcf ${reference_mag}_filtered.vcf.gz
    """
}
