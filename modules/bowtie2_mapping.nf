process BOWTIE2_MAPPING {
    conda "/opt/miniforge3/envs/mapeo_cacao"

    publishDir "${params.outdir}/bowtie2/mapping", mode: 'copy'

    debug true

    input:
    tuple val(sample_id), path(R1), path(R2), val(reference_id), val(reference_mag), path(index_files)

    output:
    tuple val(sample_id), val(reference_mag), path("${sample_id}_to_${reference_mag}.bam"), path("${sample_id}_to_${reference_mag}.bam.bai")

    script:
    """
    bowtie2 \
     --rg-id ${sample_id} \
     --rg SM:${sample_id} \
     --threads 12 \
     -x ${reference_mag} \
     -1 ${R1} \
     -2 ${R2} | \
     samtools view \
     -@12 \
     -F 4 \
     -bS - | \
     samtools sort \
     -@12 \
     -o ${sample_id}_to_${reference_mag}.bam -    
    
    samtools index -@ 12 ${sample_id}_to_${reference_mag}.bam
    """
}
