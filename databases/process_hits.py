#!/usr/bin/env python3
"""Process DIAMOND hits with variant-aware counting and proper normalization."""

import os
import sys
import pandas as pd
from collections import defaultdict

def extract_gene_id(sequence_id):
    """Extract gene ID from sequence identifier.
    Example: K00114_pwi:MWN52_00300 -> K00114"""
    return sequence_id.split('_')[0]

def count_db_variants():
    """Count number of variants per gene in the database."""
    variants_per_gene = defaultdict(int)
    
    # Read the FASTA file
    with open("kegg/target_genes.fasta") as f:
        for line in f:
            if line.startswith('>'):
                gene_id = extract_gene_id(line[1:].strip())
                variants_per_gene[gene_id] += 1
    
    return variants_per_gene

def process_diamond_hits(blast_file, total_orfs, min_identity=60, variant_identity=95):
    """Process DIAMOND hits with variant-aware counting.
    
    Args:
        blast_file: Path to DIAMOND output file (outfmt 6)
        total_orfs: Total number of ORFs in the sample
        min_identity: Minimum identity for first hit (default: 60)
        variant_identity: Minimum identity for additional variants (default: 95)
    
    Returns:
        DataFrame with normalized counts
    """
    # Read BLAST results
    cols = ['query', 'subject', 'identity', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = pd.read_csv(blast_file, sep='\t', names=cols)
    
    # Get variants per gene in database
    db_variants = count_db_variants()
    
    # Process hits
    gene_hits = defaultdict(list)
    for _, row in df.iterrows():
        gene_id = extract_gene_id(row['subject'])
        identity = row['identity']
        
        # Check if this is first hit for this gene
        if not gene_hits[gene_id]:
            if identity >= min_identity:
                gene_hits[gene_id].append(identity)
        # For additional variants, use stricter cutoff
        elif identity >= variant_identity:
            gene_hits[gene_id].append(identity)
    
    # Calculate normalized counts
    results = []
    for gene_id, hits in gene_hits.items():
        # Count number of variants found
        num_variants = len(hits)
        
        # Normalize by number of variants in database
        db_norm = num_variants / db_variants[gene_id]
        
        # Normalize by total ORFs (per million)
        final_norm = (db_norm / total_orfs) * 1_000_000
        
        results.append({
            'gene_id': gene_id,
            'variants_found': num_variants,
            'variants_in_db': db_variants[gene_id],
            'db_normalized': db_norm,
            'final_normalized': final_norm
        })
    
    return pd.DataFrame(results)

def main(blast_file, total_orfs):
    """Main function to process hits."""
    # Process hits
    results = process_diamond_hits(blast_file, total_orfs)
    
    # Sort by normalized count
    results = results.sort_values('final_normalized', ascending=False)
    
    # Save results
    output_file = blast_file.rsplit('.', 1)[0] + '_normalized.tsv'
    results.to_csv(output_file, sep='\t', index=False)
    print(f"\nResults saved to: {output_file}")
    print("\nTop 10 hits:")
    print(results.head(10))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_hits.py <blast_file> <total_orfs>")
        sys.exit(1)
    
    blast_file = sys.argv[1]
    total_orfs = int(sys.argv[2])
    main(blast_file, total_orfs)
