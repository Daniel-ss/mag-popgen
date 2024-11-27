process BOWTIE2_MAPPING {
    conda "/opt/miniforge3/envs/mapeo_cacao"

    publishDir "${params.outdir}/bowtie2/mapping", mode: 'copy'

    debug true

    input:
    tuple val(sample_id), path(R1), path(R2), val(reference_id), val(reference_mag), path(index_files)

    output:
    tuple val(sample_id), val(reference_id), path("${sample_id}_to_${reference_mag}.bam")

    script:
    """

    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`

#    echo "\$INDEX"
#    echo "${sample_id}_to_${reference_mag}.bam"     
#    touch ${sample_id}_to_${reference_mag}.bam
    bowtie2 \
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

    """
}
