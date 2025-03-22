#!/usr/bin/env python3
"""Download sequences from KEGG, COG, UniProt and NCBI."""
import os
import sys
import time
import json
import requests
from io import StringIO
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Set your email for NCBI
Entrez.email = "paulagarciamnz@gmail.com"

# Track failed downloads
FAILED_IDS_FILE = "failed_ids.json"
PROGRESS_FILE = "download_progress.json"

def load_failed_ids():
    """Load previously failed IDs."""
    if os.path.exists(FAILED_IDS_FILE):
        with open(FAILED_IDS_FILE) as f:
            return json.load(f)
    return {"kegg": [], "uniprot": [], "ncbi": []}

def save_failed_ids(failed_ids):
    """Save failed IDs for later retry."""
    with open(FAILED_IDS_FILE, 'w') as f:
        json.dump(failed_ids, f)

def load_progress():
    """Load download progress."""
    if os.path.exists(PROGRESS_FILE):
        with open(PROGRESS_FILE) as f:
            return json.load(f)
    return {"kegg": [], "uniprot": [], "ncbi": []}

def save_progress(progress):
    """Save download progress."""
    with open(PROGRESS_FILE, 'w') as f:
        json.dump(progress, f)

def get_existing_ids(fasta_file):
    """Get IDs from existing FASTA file."""
    if not os.path.exists(fasta_file):
        return set()
    
    ids = set()
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                # Extract ID from header (everything before first underscore or space)
                id = line[1:].split('_')[0].split()[0].strip()
                ids.add(id)
    return ids

def is_protein_sequence(sequence):
    """Check if a sequence is likely a protein sequence."""
    # Convert to string if it's a Seq object
    if isinstance(sequence, Seq):
        sequence = str(sequence)
    
    # Remove any whitespace and common FASTA formatting
    sequence = ''.join(sequence.split())
    
    # Check if sequence contains common amino acid letters
    amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    seq_letters = set(sequence.upper())
    
    # If more than 90% of unique characters are amino acids, it's likely a protein
    return len(seq_letters.intersection(amino_acids)) / len(seq_letters) > 0.9

def clean_sequence(seq):
    """Clean sequence by removing newlines and whitespace."""
    return ''.join(str(seq).split())

def download_kegg_sequence(kid, max_retries=3):
    """Download a sequence from KEGG with retries."""
    for attempt in range(max_retries):
        try:
            # First get the genes associated with this KO
            r = requests.get(f"https://rest.kegg.jp/link/genes/{kid}", timeout=60)
            r.raise_for_status()
            genes = [line.split()[1] for line in r.text.strip().split('\n') if line]
            
            if not genes:
                print(f"No genes found for {kid}")
                return []
            
            # Now get sequences for these genes
            sequences = []
            for gene in genes[:100]:  # Limit to 100 genes per KO
                for gene_attempt in range(max_retries):
                    try:
                        r = requests.get(f"https://rest.kegg.jp/get/{gene}/aaseq", timeout=60)
                        r.raise_for_status()
                        seq = clean_sequence(r.text.split('\n', 1)[1])
                        if is_protein_sequence(seq):
                            sequences.append(f">{kid}_{gene}\n{seq}\n")
                        break  # Success, break retry loop
                    except Exception as e:
                        if gene_attempt < max_retries - 1:
                            print(f"Error fetching {gene}: {str(e)}, retrying... ({gene_attempt + 1}/{max_retries})")
                            time.sleep(5)  # Wait longer between retries
                        else:
                            print(f"Error fetching {gene}: {str(e)}, all retries failed")
                time.sleep(0.5)  # Slightly longer delay between genes
            
            if sequences:
                print(f"Found {len(sequences)} genes", end=' ')
                return sequences
            else:
                print(f"No valid sequences found for {kid}")
                return []
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"Error: {str(e)}, retrying... ({attempt + 1}/{max_retries})")
                time.sleep(10)  # Even longer wait between full retries
            else:
                print(f"Error: {str(e)}, all retries failed")
                return []
    return []

def download_uniprot_sequence(uid, max_retries=3):
    """Download a sequence from UniProt with retries."""
    for attempt in range(max_retries):
        try:
            r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=30)
            r.raise_for_status()
            
            sequences = []
            for record in SeqIO.parse(StringIO(r.text), "fasta"):
                seq = clean_sequence(record.seq)
                if is_protein_sequence(seq):
                    sequences.append(f">{uid}\n{seq}\n")
            
            return sequences
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"UniProt error: {str(e)}, retrying... ({attempt + 1}/{max_retries})")
                time.sleep(5)
            else:
                print(f"UniProt error: {str(e)}, all retries failed")
                return []
    return []

def download_ncbi_sequence(gid, max_retries=3):
    """Download a sequence from NCBI with retries."""
    for attempt in range(max_retries):
        try:
            handle = Entrez.efetch(db="protein", id=gid, rettype="fasta", retmode="text")
            sequences = []
            for record in SeqIO.parse(handle, "fasta"):
                seq = clean_sequence(record.seq)
                if is_protein_sequence(seq):
                    sequences.append(f">{gid}\n{seq}\n")
            handle.close()
            return sequences
        except Exception as e:
            if attempt < max_retries - 1:
                print(f"NCBI error: {str(e)}, retrying... ({attempt + 1}/{max_retries})")
                time.sleep(5)
            else:
                print(f"NCBI error: {str(e)}, all retries failed")
                return []
    return []

def download_other_sequences(other_ids, output_file):
    """Download sequences from UniProt and NCBI."""
    print("\nDownloading other sequences...")
    print(f"Found {len(other_ids)} other IDs to download")
    
    # Load progress and failed IDs
    progress = load_progress()
    failed_ids = load_failed_ids()
    
    # Get existing IDs
    existing_ids = get_existing_ids(output_file)
    new_ids = [oid for oid in other_ids if oid not in existing_ids and oid not in progress["uniprot"] + progress["ncbi"]]
    print(f"Found {len(new_ids)} new IDs to download")
    
    sequences = []
    if os.path.exists(output_file):
        with open(output_file) as f:
            sequences.extend(f.readlines())
    
    for oid in new_ids:
        print(f"Fetching {oid}...", end=' ')
        # Try UniProt first
        seqs = download_uniprot_sequence(oid)
        if not seqs:
            # If not found in UniProt, try NCBI
            seqs = download_ncbi_sequence(oid)
            if seqs:
                print("Success! (NCBI)")
                progress["ncbi"].append(oid)
            else:
                print("Not found")
                failed_ids["ncbi"].append(oid)
        else:
            print("Success! (UniProt)")
            progress["uniprot"].append(oid)
        sequences.extend(seqs)
        
        # Save progress periodically
        if len(sequences) % 10 == 0:
            save_progress(progress)
            save_failed_ids(failed_ids)
            with open(output_file, 'w') as f:
                f.writelines(sequences)
    
    # Write sequences to file
    if sequences:
        with open(output_file, 'w') as f:
            f.writelines(sequences)
        print(f"Total sequences in file: {len(''.join(sequences).split('>')) - 1}")
    else:
        print("No sequences were downloaded!")

def download_kegg_sequences(kegg_ids, output_file):
    """Download sequences for KEGG IDs."""
    print("Downloading KEGG sequences...")
    print(f"Found {len(kegg_ids)} KEGG IDs to download")
    
    # Load progress and failed IDs
    progress = load_progress()
    failed_ids = load_failed_ids()
    
    # Get existing IDs
    existing_ids = get_existing_ids(output_file)
    new_ids = [kid for kid in kegg_ids if kid not in existing_ids and kid not in progress["kegg"]]
    print(f"Found {len(new_ids)} new IDs to download")
    
    sequences = []
    if os.path.exists(output_file):
        with open(output_file) as f:
            sequences.extend(f.readlines())
    
    for i, kid in enumerate(new_ids, 1):
        print(f"\nProcessing {i}/{len(new_ids)}: {kid}...", end=' ')
        seqs = download_kegg_sequence(kid)
        if seqs:
            print("Success!")
            progress["kegg"].append(kid)
        else:
            print("Failed")
            failed_ids["kegg"].append(kid)
        sequences.extend(seqs)
        
        # Save progress periodically
        if i % 5 == 0:  # Save every 5 IDs
            print(f"\nSaving progress... ({i}/{len(new_ids)} IDs processed)")
            save_progress(progress)
            save_failed_ids(failed_ids)
            with open(output_file, 'w') as f:
                f.writelines(sequences)
            print(f"Current sequences in file: {len(''.join(sequences).split('>')) - 1}")
    
    # Final save
    save_progress(progress)
    save_failed_ids(failed_ids)
    
    # Write sequences to file
    if sequences:
        with open(output_file, 'w') as f:
            f.writelines(sequences)
        print(f"\nTotal sequences in file: {len(''.join(sequences).split('>')) - 1}")
    else:
        print("\nNo sequences were downloaded!")

def retry_failed_downloads():
    """Retry downloading failed IDs."""
    print("\nRetrying failed downloads...")
    failed_ids = load_failed_ids()
    
    if not any(failed_ids.values()):
        print("No failed downloads to retry!")
        return
    
    # Retry KEGG
    if failed_ids["kegg"]:
        print(f"\nRetrying {len(failed_ids['kegg'])} failed KEGG IDs...")
        download_kegg_sequences(failed_ids["kegg"], "kegg/target_kegg.fasta")
    
    # Retry UniProt/NCBI
    other_ids = failed_ids["uniprot"] + failed_ids["ncbi"]
    if other_ids:
        print(f"\nRetrying {len(other_ids)} failed UniProt/NCBI IDs...")
        download_other_sequences(other_ids, "kegg/target_other.fasta")

def combine_and_create_diamond_db(kegg_fasta, cog_fasta, other_fasta, output_fasta):
    """Combine FASTA files and create DIAMOND database."""
    print("\nCreating combined database...")
    
    sequences = []
    
    # Add sequences from KEGG
    if os.path.exists(kegg_fasta):
        print(f"Adding sequences from {kegg_fasta}")
        with open(kegg_fasta) as f:
            sequences.extend(f.readlines())
    
    # Add sequences from COG
    if os.path.exists(cog_fasta):
        print(f"Adding sequences from {cog_fasta}")
        with open(cog_fasta) as f:
            sequences.extend(f.readlines())
    
    # Add sequences from other sources
    if os.path.exists(other_fasta):
        print(f"Adding sequences from {other_fasta}")
        with open(other_fasta) as f:
            sequences.extend(f.readlines())
    
    # Write combined sequences
    with open(output_fasta, 'w') as f:
        f.writelines(sequences)
    print(f"Total sequences in combined file: {len(''.join(sequences).split('>')) - 1}")
    
    # Create DIAMOND database
    print("\nCreating DIAMOND database...")
    os.system(f"diamond makedb --in {output_fasta} --db {output_fasta[:-6]}")

def main():
    # Get directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Input files
    kegg_ids_file = os.path.join(script_dir, 'kegg', 'target_kegg_ids.txt')
    cog_ids_file = os.path.join(script_dir, 'kegg', 'target_cog_ids.txt')
    other_ids_file = os.path.join(script_dir, 'kegg', 'target_other_ids.txt')
    
    # Output files
    kegg_fasta = os.path.join(script_dir, 'kegg', 'target_kegg.fasta')
    cog_fasta = os.path.join(script_dir, 'kegg', 'target_cog.fasta')
    other_fasta = os.path.join(script_dir, 'kegg', 'target_other.fasta')
    combined_fasta = os.path.join(script_dir, 'kegg', 'target_genes.fasta')
    
    # Read IDs
    with open(kegg_ids_file) as f:
        kegg_ids = [line.strip() for line in f if line.strip()]
    with open(other_ids_file) as f:
        other_ids = [line.strip() for line in f if line.strip()]
    
    # Download sequences
    download_kegg_sequences(kegg_ids, kegg_fasta)
    download_other_sequences(other_ids, other_fasta)
    
    # Retry any failed downloads
    retry_failed_downloads()
    
    # Combine and create DIAMOND database
    combine_and_create_diamond_db(kegg_fasta, cog_fasta, other_fasta, combined_fasta)

if __name__ == "__main__":
    main()
