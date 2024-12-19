process POGENOM {
    publishDir "${params.outdir}/pogenom", mode: 'copy'

    input:
    tuple val(sample_id), val(reference_id), path(mag), path(vcf), path(annotation)

    output:
    tuple val(sample_id), val(reference_id), path("${reference_id}"), emit: results

    script:
    """
    mkdir -p ${reference_id}

    # Download pogenom.pl script
    wget -O pogenom.pl https://raw.githubusercontent.com/EnvGen/POGENOM/refs/heads/master/pogenom.pl
    
    # Download genetic code file
    wget -O standard_genetic_code.txt https://raw.githubusercontent.com/EnvGen/POGENOM/refs/heads/master/standard_genetic_code.txt

    # Run POGENOM
    perl pogenom.pl \
     --vcf_file ${vcf} \
     --out ${reference_id}/${reference_id} \
     --gff_file ${annotation}/${sample_id}_${reference_id}.gff \
     --genetic_code_file standard_genetic_code.txt
    """
}
