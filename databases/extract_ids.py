#!/usr/bin/env python3
"""Extract KEGG, COG, and GenBank IDs from TSV file."""
import pandas as pd
import os

def extract_ids(input_file):
    """Extract KEGG, COG, and GenBank IDs from TSV file."""
    # Read TSV file
    df = pd.read_csv(input_file, sep='\t')
    
    # First get all KEGG IDs
    kegg_ids = set()
    rows_with_kegg = set()  # Keep track of rows that have KEGG IDs
    for idx, ids in df['KEGG ID'].items():
        if isinstance(ids, str):
            for id in ids.split(','):
                kid = id.strip()
                if kid and kid != 'NA' and kid.startswith('K'):
                    kegg_ids.add(kid)
                    rows_with_kegg.add(idx)
    
    # Then get COG IDs only for rows without KEGG IDs
    cog_ids = set()
    for idx, ids in df['COG ID'].items():
        if idx not in rows_with_kegg and isinstance(ids, str):
            for id in ids.split(','):
                cid = id.strip()
                if cid and cid != 'NA' and cid.startswith('COG'):
                    cog_ids.add(cid)
    
    # Extract GenBank IDs (last column 'GenBank')
    genbank_ids = set()
    if 'GenBank' in df.columns:
        for ids in df['GenBank'].dropna():
            if isinstance(ids, str) and ids.strip() != 'NA':
                for id in ids.split(','):
                    gid = id.strip()
                    if gid and gid != 'NA':
                        genbank_ids.add(gid)
    
    # Write IDs to files
    output_dir = os.path.join(os.path.dirname(input_file), 'Enceladupedia', 'databases', 'kegg')
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Found {len(kegg_ids)} KEGG IDs")
    print(f"Found {len(cog_ids)} COG IDs (only from rows without KEGG IDs)")
    print(f"Found {len(genbank_ids)} GenBank IDs")
    
    # Write KEGG IDs
    with open(os.path.join(output_dir, 'target_kegg_ids.txt'), 'w') as f:
        for kid in sorted(kegg_ids):
            f.write(f"{kid}\n")
    
    # Write COG IDs
    with open(os.path.join(output_dir, 'target_cog_ids.txt'), 'w') as f:
        for cid in sorted(cog_ids):
            f.write(f"{cid}\n")
    
    # Write GenBank IDs
    with open(os.path.join(output_dir, 'target_other_ids.txt'), 'w') as f:
        for gid in sorted(genbank_ids):
            f.write(f"{gid}\n")

if __name__ == "__main__":
    input_file = "/Users/paulagarciamartinez/Desktop/KAUST_2_SEM/Master_Thesis/RS_HV_Metagenomics/target_genes/All_target_genes_Icy_Moons_TSV_cleaned.txt"
    extract_ids(input_file)
