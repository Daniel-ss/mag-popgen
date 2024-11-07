process EXTRACT_CONTIG_NAMES {
    conda "/opt/miniforge3/envs/AI"

    publishDir "${params.outdir}/contig_names", mode: 'copy'
    
    input:
    tuple val(sample_id), path(mag)
    
    output:
    tuple val(sample_id), path("${sample_id}_${mag.baseName}_contigs.txt"), emit: contig_names
    
    script:
    """
    python ${projectDir}/bin/process_mags.py ${mag} ${sample_id}_${mag.baseName}_contigs.txt
    """
}
