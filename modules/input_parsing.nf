process PARSE_INPUTS {
    input:
    path bam_paths
    path mag_paths
    
    output:
    path "parsed_inputs.json"
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import sys
    
    def parse_file(file_path):
        data = {}
        try:
            with open(file_path, 'r') as f:
                header = f.readline().strip().split('\t')
                for line in f:
                    fields = line.strip().split('\t')
                    row_data = dict(zip(header, fields))
                    sample_id = row_data['sample_id']
                    if sample_id not in data:
                        data[sample_id] = []
                    data[sample_id].append(row_data)
        except Exception as e:
            print(f"Error parsing file {file_path}: {str(e)}", file=sys.stderr)
            sys.exit(1)
        return data
    
    try:
        bam_data = parse_file("${bam_paths}")
        mag_data = parse_file("${mag_paths}")
        
        # Verify data structure
        for sample_id in mag_data:
            for mag in mag_data[sample_id]:
                assert 'mag_id' in mag and 'mag_path' in mag, f"Missing required fields in MAG data for sample {sample_id}"
        
        for sample_id in bam_data:
            for bam in bam_data[sample_id]:
                assert 'bam_path' in bam, f"Missing bam_path in BAM data for sample {sample_id}"
        
        parsed_data = {
            'bam_files': bam_data,
            'mag_files': mag_data
        }
        
        with open('parsed_inputs.json', 'w') as f:
            json.dump(parsed_data, f, indent=2)
            
    except Exception as e:
        print(f"Error processing input files: {str(e)}", file=sys.stderr)
        sys.exit(1)
    """
}
