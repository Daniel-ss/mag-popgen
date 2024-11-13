process BOWTIE2_MAPPING {
    conda "/opt/miniforge3/envs/mapeo_cacao"

    publishDir "${params.outdir}/bowtie2/mapping", mode: 'copy'

    input:
    tuple val(sample_id), path(mag)
    path(index)

    output:
    tuple val(sample_id), path("${sample_id}_${mag.baseName}*"), emit: mag_bam_file

    script:
    """
    echo "${index}"
    echo "${sample_id}_${mag.baseName}"
#    bowtie2 \
#     --threads 12 \
#     -x test_data/COMEBin_106 \
#     -1 ~/Documents/Tesis/Data/${sample}_1_clean.fastq.gz \
#     -2 ~/Documents/Tesis/Data/${sample}_2_clean.fastq.gz | \
#     samtools view \
#     -@12 \
#     -F 4 \
#     -bS - | \
#     samtools sort \
#     -@12 \
#     -o test_data/${sample}_to_COMEBIN_106.bam -    

    """
}
