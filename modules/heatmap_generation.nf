process GENERATE_HEATMAP {
    conda "/opt/miniforge3/envs/AI"
    publishDir "${params.outdir}/heatmaps", mode: 'copy'
    
    input:
    tuple val(sample_id), path(coverage_file)
    
    output:
    path "${sample_id}_abundance_heatmap.png"
    
    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Read the coverage file
    Abundance = pd.read_csv("${coverage_file}", sep="\t")
    
    # Set the 'Genome' column as the index
    Abundance.set_index('Genome', inplace=True)
    
    # Remove unmapped reads
    data = Abundance[Abundance.index != 'unmapped']
    
    # Identify sample columns (those containing 'Relative Abundance')
    sample_columns = [col for col in data.columns if 'Relative Abundance' in col]
    
    # Select only the sample columns
    abundance_data = data[sample_columns]
    
    # Clean up column names to just sample names
    abundance_data.columns = [col.split('_to_')[0] for col in abundance_data.columns]
    
    # Create a mask for unmapped rows
    unmapped_mask = np.zeros_like(abundance_data, dtype=bool)
    unmapped_mask[abundance_data.index == 'unmapped'] = True
    
    # Create a copy of the data for coloration
    heatmap_data = abundance_data.copy()
    
    # Replace unmapped values with NaN for coloration
    heatmap_data.loc[unmapped_mask] = np.nan
    
    #Figure size and label size
    figure_width = len(abundance_data.columns)*2
    figure_height = len(abundance_data)*0.5
    label_fontsize = min(figure_width, figure_height) * 1.5   

    # Create the figure and axes
    plt.figure(figsize=(len(abundance_data.columns)*2, len(abundance_data)*0.5))
    
    # Create the heatmap with a custom color handling
    ax = sns.heatmap(heatmap_data, 
                     annot=abundance_data,  # Use original data for annotations
                     fmt='.2f',  # Format to show 2 decimal places
                     cmap='YlOrRd',  # Color palette 
                     cbar_kws={'label': 'Relative Abundance (%)'},
                     linewidths=0.5,  # Add lines between cells
                     mask=unmapped_mask,  # Mask unmapped rows from coloration
                     annot_kws={"color": "black"})  # Ensure annotations are always black
    
    ax.figure.axes[-1].set_ylabel('Relative abundance (%)', size=label_fontsize)

    plt.xlabel('Samples')
    plt.ylabel('Genomes')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Save the plot
    plt.savefig('${sample_id}_abundance_heatmap.png', dpi=300, bbox_inches='tight')
    """
}
