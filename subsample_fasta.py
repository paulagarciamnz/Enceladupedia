#!/usr/bin/env python3

from Bio import SeqIO
import sys
import math

def subsample_fasta(input_file, output_file, percentage=5):
    """
    Take a percentage of sequences from a FASTA file while preserving order.
    
    Args:
        input_file (str): Path to input FASTA file
        output_file (str): Path to output FASTA file
        percentage (float): Percentage of sequences to keep (default: 5)
    """
    # Count total sequences
    total_sequences = sum(1 for _ in SeqIO.parse(input_file, "fasta"))
    num_sequences_to_keep = math.ceil(total_sequences * (percentage / 100))
    
    print(f"Total sequences: {total_sequences}")
    print(f"Keeping {num_sequences_to_keep} sequences ({percentage}%)")
    
    # Calculate step size to evenly sample throughout the file
    step = total_sequences / num_sequences_to_keep
    
    # Read and write sequences
    sequences_to_write = []
    for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
        if i % step < 1:  # This ensures we take evenly spaced sequences
            sequences_to_write.append(record)
    
    # Write selected sequences
    SeqIO.write(sequences_to_write, output_file, "fasta")
    print(f"Written {len(sequences_to_write)} sequences to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python subsample_fasta.py <input_fasta>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = input_file.rsplit('.', 1)[0] + "_5percent.fna"
    subsample_fasta(input_file, output_file)
