include { CHECKM2 } from '../modules/checkm2'
include { PARSE_INPUTS} from '../modules/input_parsing'
include { COVERM } from '../modules/coverm'
include { EXTRACT_CONTIG_NAMES } from '../modules/mag_to_contig'
include { GENERATE_HEATMAP } from '../modules/heatmap_generation'
include { BOWTIE2_BUILD } from '../modules/bowtie2_build'
include { BOWTIE2_MAPPING } from '../modules/bowtie2_mapping'
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
    //EXTRACT_CONTIG_NAMES.out.contig_names.view { "Contig names output: $it" }    

    // Group BAM files by reference genome
    grouped_bams = bam_ch
    .map { sample_id, reference_id, bam -> tuple(reference_id, bam) }
    .groupTuple(by: 0)

    // Combine MAGs with grouped BAM files
    mag_bam_combined = mag_ch
    .map { sample_id, mag -> tuple(mag.baseName, mag) }
    .join(grouped_bams)

    // Run COVERM process
 
    // Parse input files
    PARSE_INPUTS(params.bam_paths, params.mag_paths)

    // Read and process the parsed inputs
    parsed_data = PARSE_INPUTS.out
        .map { file -> 
            def json_text = file.text
            def slurper = new groovy.json.JsonSlurper()
            return slurper.parseText(json_text)
        }

    // Create a channel for each sample
    sample_ch = parsed_data.flatMap { data ->
        println "Processing parsed data..."
        
        data.mag_files.collect { sample_id, mags ->
            println "Processing sample: ${sample_id}"
            //println "MAGs: ${mags}"
            
            try {
                // Create genome definition file
                def genome_def = "genome_def_${sample_id}.txt"
                def genome_def_content = new StringBuilder()
                
                // Process MAGs
                mags.each { mag ->
                    //println "Processing MAG: ${mag}"
                    if (!mag.mag_path || !mag.mag_id) {
                        println "Warning: Missing mag_path or mag_id in ${mag}"
                        return
                    }
                    
                    def mag_file = file(mag.mag_path)
                    if (mag_file.exists()) {
                        mag_file.eachLine { line ->
                            if (line.startsWith('>')) {
                                def contig_name = line.substring(1).split()[0]
                                genome_def_content.append("${mag.mag_id}\t${contig_name}\n")
                            }
                        }
                    } else {
                        println "Warning: MAG file not found: ${mag.mag_path}"
                    }
                }
                
                // Get corresponding BAM files
                def bams = data.bam_files[sample_id]
                if (!bams) {
                    println "Warning: No BAM files found for sample ${sample_id}"
                    return null
                }
                
                // Collect BAM paths
                def bam_paths = bams.collect { bam -> 
                    if (!bam.bam_path) {
                        println "Warning: null BAM path in ${bam}"
                        return null
                    }
                    return file(bam.bam_path)
                }.findAll { it != null }
                
                if (bam_paths.isEmpty()) {
                    println "Warning: No valid BAM paths for sample ${sample_id}"
                    return null
                }
                
                // Only create output if we have content
                if (genome_def_content.length() > 0) {
                    def genome_def_file = file("${workflow.workDir}/genome_def_${sample_id}.txt")
                    genome_def_file.text = genome_def_content.toString()
                    
                    //println "Created genome definition file for ${sample_id}"
                    //println "BAM paths: ${bam_paths}"
                    
                    return tuple(sample_id, genome_def_file, bam_paths)
                }
                
                return null
            } catch (Exception e) {
                println "Error processing sample ${sample_id}: ${e.message}"
                e.printStackTrace()
                return null
            }
        }.findAll { it != null }
    }
    
    // Add verification before COVERM
    //sample_ch.subscribe { println "Sample channel entry: $it" }
 
    COVERM(sample_ch)

    // Generate heatmaps for each coverage file
    GENERATE_HEATMAP(COVERM.out.coverage)

    BOWTIE2_BUILD(mag_ch)

    BOWTIE2_MAPPING(mag_ch, BOWTIE2_BUILD.out.flatten())

    FREEBAYES(mag_bam_combined)

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
