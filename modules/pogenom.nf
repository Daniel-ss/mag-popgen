process POGENOM {
    
    publishDir "${params.outdir}/pogenom", mode: 'copy'

    input:
    tuple val(sample_id), val(reference_id), path(mag), path(bam), path(vcf), path(coverage)
    path(gff_file)

    output:
    tuple val(sample_id), val(reference_id), path("${sample_id}_${reference_id}_pogenom_output.txt"), emit: results

    script:
    """
    # Download pogenom.pl script
    wget -O pogenom.pl https://raw.githubusercontent.com/EnvGen/POGENOM/refs/heads/master/pogenom.pl
    
    # Download genetic code file
    wget -O standard_genetic_code.txt https://raw.githubusercontent.com/EnvGen/POGENOM/refs/heads/master/standard_genetic_code.txt
    
    # Run POGENOM
    perl pogenom.pl \
     --vcf_file ${vcf} \
     --out ${sample_id}_${reference_id}_pogenom_output.txt \
     --gff_file ${gff_file} \
     --genetic_code_file standard_genetic_code.txt
    """
}
