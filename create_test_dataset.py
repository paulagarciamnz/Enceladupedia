from Bio import SeqIO
import random
import os

def create_subset(input_file, output_file, percentage=5):
    # Read all records
    records = list(SeqIO.parse(input_file, "fasta"))
    
    # Calculate number of sequences to keep
    num_sequences = len(records)
    num_to_keep = int(num_sequences * (percentage / 100))
    
    # Randomly select sequences
    selected_records = random.sample(records, num_to_keep)
    
    # Write selected sequences to new file
    SeqIO.write(selected_records, output_file, "fasta")
    print(f"Created subset with {num_to_keep} sequences out of {num_sequences} total sequences")

if __name__ == "__main__":
    # Get the parent directory path
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Define input and output paths
    input_file = os.path.join(parent_dir, "3300075498_contigs.fna")
    output_file = "test_subset.fna"
    
    if not os.path.exists(input_file):
        print(f"Error: Input file not found at {input_file}")
        exit(1)
        
    create_subset(input_file, output_file)
