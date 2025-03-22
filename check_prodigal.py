#!/usr/bin/env python3

import subprocess
import os
import re

# Path to the test file
test_file = "uploads/test_subset.fna"

# Get only the first 1000 sequences for testing
output_first_1000 = "uploads/first_1000.fna"

def count_fasta_sequences(file_path):
    count = 0
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                count += 1
    return count

# Create a file with only the first 1000 sequences
with open(test_file, 'r') as f_in, open(output_first_1000, 'w') as f_out:
    seq_count = 0
    writing = True
    
    for line in f_in:
        if line.startswith('>'):
            seq_count += 1
            if seq_count > 1000:
                writing = False
        
        if writing:
            f_out.write(line)

print(f"Created file with the first 1000 sequences out of {count_fasta_sequences(test_file)}")

# Run Prodigal on the first 1000 sequences
prodigal_out = "uploads/prodigal_test.faa"
prodigal_cmd = f"prodigal -i {output_first_1000} -a {prodigal_out} -p meta -q"

print(f"Running command: {prodigal_cmd}")
process = subprocess.Popen(prodigal_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
stdout, stderr = process.communicate()

if process.returncode != 0:
    print(f"Error running Prodigal: {stderr}")
else:
    print("Prodigal finished successfully")
    
    # Count the number of sequences in the output file
    protein_count = 0
    with open(prodigal_out, 'r') as f:
        for line in f:
            if line.startswith('>'):
                protein_count += 1
    
    print(f"Output file contains {protein_count} protein sequences")
