process BOWTIE2_MAPPING {
    conda "/opt/miniforge3/envs/mapeo_cacao"

    publishDir "${params.outdir}/bowtie2/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(mag), path(R1), path(R2) 

    output:
    tuple val(sample_id), path("${sample_id}_${R1}_${R2}_${mag.baseName}*"), emit: mag_bam_file

    script:
    """
#    touch ${sample_id}_${R1}_${R2}_${mag.baseName}.bam
#    touch ${sample_id}_to_${mag.baseName}.txt


    bowtie2 \
     --threads 12 \
     -x ${mag.baseName} \
     -1 ${R1} \
     -2 ${R2} | \
     samtools view \
     -@12 \
     -F 4 \
     -bS - | \
     samtools sort \
     -@12 \
     -o ${sample_id}_to_${mag.baseName}.bam -    

    """
}
