
<img width="970" alt="image" src="https://github.com/user-attachments/assets/08d6fc47-8ed0-490e-9e3c-ae6df4298441" />

## Quick assess the astrobiology potential of your metagenomic data
# Enceladupedia 🪐🦠🧫
A web application for analyzing metagenomic sequences and identifying target genes related to methanogenesis and other metabolic pathways.


## Getting Started

Enceladupedia works through a user-friendly web interface, but you'll need to set up the environment using command-line instructions first. If you're not familiar with command lines, don't worry! The process is straightforward.

We all have started somewhere, you might find this resources useful:
- [Introduction to Command Line](https://tutorial.djangogirls.org/en/intro_to_command_line/)
- [Installing Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/Enceladupedia.git
cd Enceladupedia
```

2. Create and activate a conda environment:
```bash
conda create -n enceladupedia python=3.9
conda activate enceladupedia
```

3. Install dependencies:
```bash
pip install -r requirements.txt
conda install -c bioconda diamond prodigal
```

4. Download and prepare the KEGG database:
```bash
# Create database directory
mkdir -p databases/kegg
cd databases/kegg

# Download UniProt data with KEGG mappings
# Note: The large database file (idmapping_selected.tab.gz) is not included in the repository due to size constraints
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz
curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# Extract files
gunzip idmapping_selected.tab.gz
gunzip uniprot_sprot.fasta.gz

# Process files to create KEGG-specific FASTA
python prepare_kegg_db.py
```

5. Build the custom DIAMOND database with GenBank entries:
```bash
# Return to the main directory
cd ../../

python databases/download_sequences.py

# The download_sequences.py script will create a DIAMOND database including all our target genes

```

6. Run the application:
```bash
python app.py
```

The application will be available at http://localhost:8888

## Features

- Upload and analyze assembled sequences (optimized for contigs, you might want to change parameters in the app.py to work with genomes)
- Identify genes using Prodigal (optimized for prokaryotic metagenomes)
- Match predicted proteins against KEGG database using DIAMOND
- Generate interactive visualizations of astrobiology-relevant features
- Cache results for faster repeated analyses

## Requirements

- Python 3.9+
- DIAMOND aligner
- Prodigal gene finder
- ~20GB disk space for databases
- 4GB+ RAM recommended

## Database Setup Details

The application uses a local DIAMOND database created from UniProt-KEGG mappings. The setup process:

1. Downloads UniProt's SwissProt database (~2GB compressed)
2. Downloads UniProt-KEGG ID mappings (~1GB compressed)
3. Processes these to create a KEGG-specific FASTA file
4. Creates a DIAMOND database for fast searching

This process:
- Takes 1-2 hours depending on your internet connection
- Requires ~20GB temporary disk space
- Results in a ~5GB DIAMOND database

## Alternative Database Options

If you can't download the full database, you can:

1. Use a smaller subset of KEGG (instructions in `docs/small_db_setup.md`)
2. Use pre-computed DIAMOND database (contact maintainers)
3. Use the online version at [your-website-url]

## Contributing

Contributions are welcome! Please read `CONTRIBUTING.md` for guidelines.

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.

## Citation

If you use this tool in your research, please cite:
[Your citation information]
