#!/usr/bin/env python3
"""Clean up TSV file by fixing formatting issues."""
import pandas as pd
import numpy as np

def clean_tsv(input_file, output_file):
    """Clean up TSV file by fixing formatting issues."""
    # Read the file as raw text first to count fields
    with open(input_file, 'r') as f:
        header = f.readline().strip().split('\t')
        expected_fields = len(header)
        print(f"Expected {expected_fields} fields based on header")
    
    # Read file with a more lenient parser
    df = pd.read_csv(input_file, sep='\t', header=0, dtype=str, on_bad_lines='warn', quoting=3)
    print(f"Read {len(df)} rows")
    
    # Clean up each column
    for col in df.columns:
        # Replace NaN with 'NA'
        df[col] = df[col].fillna('NA')
        
        # Remove any newlines and extra spaces
        df[col] = df[col].apply(lambda x: ' '.join(str(x).split()))
        
        # Remove any quotes
        df[col] = df[col].str.replace('"', '')
        
        # Strip whitespace
        df[col] = df[col].str.strip()
    
    # Save cleaned file
    df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    print(f"Saved cleaned file to {output_file}")

if __name__ == "__main__":
    input_file = "/Users/paulagarciamartinez/Desktop/KAUST_2_SEM/Master_Thesis/RS_HV_Metagenomics/target_genes/All_target_genes_Icy_Moons_TSV.txt"
    output_file = "/Users/paulagarciamartinez/Desktop/KAUST_2_SEM/Master_Thesis/RS_HV_Metagenomics/target_genes/All_target_genes_Icy_Moons_TSV_cleaned.txt"
    clean_tsv(input_file, output_file)
