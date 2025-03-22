#!/usr/bin/env python3
"""
Process UniProt files to create a KEGG-specific FASTA file.
This script:
1. Reads the idmapping file to get UniProt to KEGG mappings
2. Extracts relevant sequences from Swiss-Prot
3. Creates a FASTA file with only KEGG-mapped sequences
"""

import sys
from collections import defaultdict

def read_kegg_mappings(mapping_file):
    """Read UniProt to KEGG mappings."""
    print("Reading KEGG mappings...")
    kegg_map = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                uniprot_id = parts[0]
                # KEGG IDs are in column 7
                if len(parts) > 7 and parts[7]:  # Check if KEGG ID exists
                    kegg_ids = parts[7].split('; ')
                    kegg_map[uniprot_id] = kegg_ids
    print(f"Found {len(kegg_map)} UniProt IDs with KEGG mappings")
    return kegg_map

def process_fasta(fasta_file, kegg_map, output_file):
    """Extract sequences with KEGG mappings and write to new FASTA."""
    print("Processing FASTA file...")
    current_id = None
    current_seq = []
    sequences_written = 0
    
    with open(output_file, 'w') as out:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Write previous sequence if it had KEGG mapping
                    if current_id and current_id in kegg_map:
                        for kegg_id in kegg_map[current_id]:
                            out.write(f">{kegg_id}|{current_id}\n")
                            out.write("".join(current_seq))
                            sequences_written += 1
                    
                    # Start new sequence
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Don't forget the last sequence
            if current_id and current_id in kegg_map:
                for kegg_id in kegg_map[current_id]:
                    out.write(f">{kegg_id}|{current_id}\n")
                    out.write("".join(current_seq))
                    sequences_written += 1
    
    print(f"Wrote {sequences_written} sequences to {output_file}")

def main():
    # Read KEGG mappings
    kegg_map = read_kegg_mappings("idmapping_selected.tab")
    
    # Process FASTA file
    process_fasta("uniprot_sprot.fasta", kegg_map, "kegg_proteins.fasta")
    
    print("Done! Now you can create a DIAMOND database with:")
    print("diamond makedb --in kegg_proteins.fasta --db kegg_proteins")

if __name__ == "__main__":
    main()
