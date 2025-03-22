#!/usr/bin/env python3
"""Get descriptions for KEGG modules and combine with custom modules."""
import os
import requests
import pandas as pd
import time

def get_kegg_module_description(mid):
    """Get description for a KEGG module ID."""
    try:
        # Be nice to the KEGG API
        time.sleep(1)
        
        # Get module info
        r = requests.get(f"https://rest.kegg.jp/get/{mid}", timeout=10)
        r.raise_for_status()
        
        # Parse the response
        name = ""
        definition = ""
        lines = r.text.split('\n')
        for line in lines:
            if line.startswith('NAME'):
                name = line.split('NAME', 1)[1].strip()
            elif line.startswith('DEFINITION'):
                definition = line.split('DEFINITION', 1)[1].strip()
                break
        
        # Return name and definition if available
        if name and definition:
            return f"{name} - {definition}"
        elif name:
            return name
        elif definition:
            return definition
        return "No description available"
    except Exception as e:
        print(f"Error fetching {mid}: {str(e)}")
        return "Error fetching description"

def main():
    # Read the modules file
    modules_file = "/Users/paulagarciamartinez/Desktop/KAUST_2_SEM/Master_Thesis/RS_HV_Metagenomics/target_genes/modules.txt"
    with open(modules_file) as f:
        modules = [line.strip() for line in f if line.strip()]
    
    # Read the target genes file to get custom module descriptions
    genes_file = "/Users/paulagarciamartinez/Desktop/KAUST_2_SEM/Master_Thesis/RS_HV_Metagenomics/target_genes/All_target_genes_Icy_Moons_TSV.txt"
    genes_df = pd.read_csv(genes_file, sep='\t')
    
    # Create a dictionary to store module descriptions
    descriptions = {}
    
    # Process each module
    for module in modules:
        if module.startswith('M'):
            # Get KEGG module description
            desc = get_kegg_module_description(module)
            descriptions[module] = desc
        else:
            # For custom modules, look up in the metabolism column
            matching_rows = genes_df[genes_df['metabolism'].str.lower() == module.lower().replace('_', ' ')]
            if not matching_rows.empty:
                # Get unique pathways
                pathways = matching_rows['pathway'].unique()
                # Join multiple pathways if they exist
                custom_desc = ' | '.join([p for p in pathways if isinstance(p, str)])
                if not custom_desc:
                    custom_desc = module.replace('_', ' ')
            else:
                custom_desc = module.replace('_', ' ')
            descriptions[module] = custom_desc
    
    # Create output DataFrame
    df = pd.DataFrame(list(descriptions.items()), columns=['Module', 'Description'])
    
    # Save to TSV
    output_file = os.path.join(os.path.dirname(modules_file), 'module_descriptions.tsv')
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Saved module descriptions to: {output_file}")
    print("\nModule descriptions:")
    print(df.to_string(index=False))

if __name__ == "__main__":
    main()
