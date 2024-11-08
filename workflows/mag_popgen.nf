include { CHECKM2 } from '../modules/checkm2'
include { COVERM } from '../modules/coverm'
include { EXTRACT_CONTIG_NAMES } from '../modules/mag_to_contig'
include { FREEBAYES } from '../modules/freebayes'
include { PROKKA } from '../modules/prokka'
include { POGENOM } from '../modules/pogenom'

workflow MAG_POPGEN {
    // Create channels for input data
    mag_ch = Channel
        .fromPath(params.mag_paths)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample_id, file(row.mag_path)) }

    bam_ch = Channel
        .fromPath(params.bam_paths)
        .splitCsv(header:true, sep:'\t')
        .map { row -> tuple(row.sample_id, row.reference_id, file(row.bam_path)) }

    metadata_ch = Channel.fromPath(params.metadata)
    
    // Combine MAG and BAM channels
    combined_ch = mag_ch.combine(bam_ch, by: 0)
        .map { sample_id, mag, reference_id, bam -> 
            reference_id == mag.baseName ? tuple(sample_id, reference_id, mag, bam) : null 
        }
        .filter { it != null }

    // Run processes
    CHECKM2(mag_ch)
    EXTRACT_CONTIG_NAMES(mag_ch)
    EXTRACT_CONTIG_NAMES.out.contig_names.view { "Contig names output: $it" }    

    // Group BAM files by reference genome
    grouped_bams = bam_ch
    .map { sample_id, reference_id, bam -> tuple(reference_id, bam) }
    .groupTuple(by: 0)

    // Combine MAGs with grouped BAM files
    mag_bam_combined = mag_ch
    .map { sample_id, mag -> tuple(mag.baseName, mag) }
    .join(grouped_bams)

    FREEBAYES(mag_bam_combined)
    
    // Combine MAGs with grouped BAM files and contig names
    //coverm_input = mag_ch
    //    .map { sample_id, mag -> tuple(mag.baseName, mag, sample_id) }
    //   .combine(grouped_bams, by: 0)
    //    .combine(EXTRACT_CONTIG_NAMES.out.contig_names, by: 0)
    //    .map { mag_basename, mag, sample_id, bams, contig_names -> 
    //        tuple(mag_basename, mag, sample_id, bams.flatten(), contig_names)
    //    }

    coverm_input = mag_bam_combined
     .join(EXTRACT_CONTIG_NAMES.out.contig_names, by: 0)
    

    // Run COVERM process
    log.info "About to execute COVERM process"
    coverm_input.view { "### COVERM input ### $it" }
    COVERM(coverm_input)

    PROKKA(mag_ch)
    
    pogenom_input = combined_ch
        .join(FREEBAYES.out.vcf, by: [0, 1])
        .join(COVERM.out.coverage, by: [0, 1])
    POGENOM(pogenom_input)

    //Collect and emit results
    results = CHECKM2.out.results
       .mix(COVERM.out.coverage)
       .mix(FREEBAYES.out.vcf)
       .mix(PROKKA.out.annotation)
       .mix(POGENOM.out.results)
       .collect()

    emit:
    results
}
