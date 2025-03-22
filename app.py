import os
import sys
import time
import logging
import traceback
import subprocess
import tempfile
import pandas as pd
from threading import Thread
from flask import Flask, request, jsonify, render_template, make_response, g
from werkzeug.utils import secure_filename
from Bio import SeqIO
import csv
from io import StringIO
import re
import random
import json

# Path to module descriptions file - made absolute
module_file = '/Users/paulagarciamartinez/Desktop/KAUST_2_SEM/Master_Thesis/RS_HV_Metagenomics/target_genes/module_descriptions.txt'  # Now reading from .txt with Enceladus and Europa scores

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Initialize app
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev_key_for_flask')
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uploads')
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB limit

# Initialize persistent storage for results
try:
    if not os.path.exists(app.config['UPLOAD_FOLDER']):
        os.makedirs(app.config['UPLOAD_FOLDER'])
except Exception as e:
    logging.error(f"Error creating upload folder: {e}")

@app.before_request
def before_request():
    """Initialize request globals."""
    # Initialize results storage if not already set
    if not hasattr(g, 'analysis_results'):
        g.analysis_results = None
    
    # Try to load saved results only if analysis_results is not set
    if g.analysis_results is None:
        results_path = os.path.join(app.config['UPLOAD_FOLDER'], 'module_results.json')
        if os.path.exists(results_path):
            try:
                with open(results_path, 'r') as f:
                    g.analysis_results = json.load(f)
            except Exception as e:
                logging.error(f"Error loading saved results: {e}")

# Configure Flask to not cache templates
app.config['TEMPLATES_AUTO_RELOAD'] = True
app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0

# Global progress tracking
progress = {
    'status': 'idle',
    'stage': None,
    'stage_text': 'Waiting for file...',
    'stage_progress': 0,
    'gene_count': 0,
    'error': None,
    'results': None,
    'start_time': None,
    'current_step_start_time': None,
    'total_sequences': 0,
    'processed_sequences': 0,
    'done': False
}

gene_module_mapping = {}

def run_prodigal(input_file):
    """Run Prodigal gene prediction."""
    try:
        # Create temporary files for output
        with tempfile.NamedTemporaryFile(suffix='.faa', delete=False) as proteins_file:
            # Run Prodigal
            cmd = [
                'prodigal',
                '-i', input_file,
                '-a', proteins_file.name,  # Protein translations
                '-p', 'meta',             # Metagenomic mode
                '-m',                     # Treat runs of N as masked sequence
                '-n',                     # Bypass Shine-Dalgarno trainer
                '-q'                      # Run quietly
            ]
            
            logging.info(f"Running Prodigal command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise Exception(f"Prodigal failed: {result.stderr}")
            
            # Read protein sequences
            protein_seqs = []
            for record in SeqIO.parse(proteins_file.name, 'fasta'):
                protein_seqs.append(str(record.seq).rstrip('*'))
            
            gene_count = len(protein_seqs)
            logging.info(f"Prodigal found {gene_count} genes")
            return protein_seqs, gene_count
            
    except Exception as e:
        logging.error(f"Error running Prodigal: {str(e)}")
        raise
    finally:
        # Cleanup
        if 'proteins_file' in locals():
            os.unlink(proteins_file.name)

def batch_diamond_search(protein_seqs, threads=4):
    """Run DIAMOND search on the protein sequences."""
    logging.info(f"Running DIAMOND search with {len(protein_seqs)} sequences using {threads} threads")
    
    # Create a temporary file for the protein sequences
    with tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False) as temp_faa:
        for i, seq in enumerate(protein_seqs):
            temp_faa.write(f">gene_{i}\n{seq}\n")
        temp_faa_path = temp_faa.name
    
    # Create a temporary file for the DIAMOND output
    temp_tsv = tempfile.NamedTemporaryFile(suffix='.tsv', delete=False)
    temp_tsv_path = temp_tsv.name
    temp_tsv.close()
    
    try:
        # Define DIAMOND command
        cmd = [
            'diamond', 'blastp',
            '--db', 'databases/kegg/target_genes',
            '--query', temp_faa_path,
            '--out', temp_tsv_path,
            '--outfmt', '6',
            '--max-target-seqs', '100',
            '--threads', str(threads),
            '--more-sensitive',
            '--id', '50'
        ]
        
        # Run DIAMOND
        logging.info(f"Running DIAMOND command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error(f"DIAMOND search failed: {result.stderr}")
            raise Exception(f"DIAMOND search failed: {result.stderr}")
        
        # Read DIAMOND output
        cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        df = pd.read_csv(temp_tsv_path, sep='\t', names=cols)
        
        logging.info(f"DIAMOND search completed with {len(df)} hits")
        return df
    
    except Exception as e:
        logging.error(f"Error in DIAMOND search: {str(e)}")
        # Return empty DataFrame in case of error - we could raise instead
        return pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                                    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
    
    finally:
        # Clean up temporary files
        if os.path.exists(temp_faa_path):
            os.unlink(temp_faa_path)
        if os.path.exists(temp_tsv_path):
            os.unlink(temp_tsv_path)

def process_sequence_results(hits_df):
    """Process the sequence results from DIAMOND to extract matched genes.
    
    Args:
        hits_df: DataFrame with DIAMOND hits
        
    Returns:
        Tuple containing: (list of gene IDs found in hits, count of genes with functions, mapping of gene IDs to query IDs)
    """
    logging.info(f"Processing results with {len(hits_df)} hits")
    
    # Extract the gene IDs from the subject identifiers
    # Handle different formats: KO IDs (K00114), GenBank accessions, and various separator formats
    gene_hits = set()
    gene_ids_to_query_ids = {}
    
    for _, row in hits_df.iterrows():
        subject_id = row['sseqid'] if 'sseqid' in row else row.get('subject_id', '')
        query_id = row['qseqid'] if 'qseqid' in row else row.get('query_id', '')
        gene_id = None
        
        # Case 1: Format KXXXXX_organism:gene
        if '_' in subject_id and subject_id.startswith('K'):
            ko_id = subject_id.split('_')[0]
            if len(ko_id) >= 6:  # KO IDs like K00001
                gene_id = ko_id
        
        # Case 2: Format with | separators (old format)
        elif '|' in subject_id:
            parts = subject_id.split('|')
            if parts and len(parts) > 1:
                potential_gene = parts[-1]
                if potential_gene.startswith('K') and len(potential_gene) >= 6:
                    gene_id = potential_gene
                else:
                    # This might be a GenBank ID
                    gene_id = potential_gene
        
        # Case 3: Direct KO IDs without separators
        elif subject_id.startswith('K') and len(subject_id) >= 6 and not '_' in subject_id:
            gene_id = subject_id
        
        # Case 4: GenBank ID or other format (specific patterns we're looking for)
        elif (re.match(r'^[A-Z0-9]+\.[0-9]+$', subject_id) or  # Accession.version (e.g., ABG91829.1)
              re.match(r'^[A-Z][A-Z][A-Z][0-9][0-9]_[0-9]+$', subject_id) or  # Locus tag (e.g., EGC82_06255)
              re.match(r'^[A-Z][0-9][A-Z][0-9][0-9]_[0-9]+$', subject_id) or  # Locus tag variation
              re.match(r'^[A-Z][0-9][A-Z][A-Z][0-9]_[0-9]+$', subject_id)):   # Another locus tag variation
            gene_id = subject_id
        
        # Case 5: Any other ID - assume it could be a valid ID
        else:
            gene_id = subject_id
        
        # Store the gene ID if found
        if gene_id:
            gene_hits.add(gene_id)
            if gene_id not in gene_ids_to_query_ids:
                gene_ids_to_query_ids[gene_id] = []
            if query_id not in gene_ids_to_query_ids[gene_id]:
                gene_ids_to_query_ids[gene_id].append(query_id)
    
    # Count unique genes that have functions
    genes_with_functions = set()
    for gene_id_list in gene_ids_to_query_ids.values():
        genes_with_functions.update(gene_id_list)
    
    logging.info(f"Found {len(gene_hits)} unique gene identifiers")
    logging.info(f"Found {len(genes_with_functions)} genes with functions")
    
    return list(gene_hits), len(genes_with_functions), gene_ids_to_query_ids

def extract_ko_ids(definition):
    """Extract KO IDs from a KEGG module definition.
    
    KEGG module definitions use a complex notation with +, -, and parentheses.
    This function extracts all KO IDs regardless of their relationship.
    
    Args:
        definition: String containing the KEGG module definition
    
    Returns:
        List of KO IDs found in the definition
    """
    if not definition:
        return []
    
    # Use regex to extract all KO IDs (K followed by 5 digits)
    ko_pattern = r'K\d{5}'
    ko_ids = re.findall(ko_pattern, definition)
    
    # Return unique KO IDs
    return list(set(ko_ids))

def extract_kegg_module_genes(definition):
    """Extract KO IDs from a KEGG module definition.
    
    KEGG module definitions use a complex notation with +, -, and parentheses.
    This function extracts all KO IDs regardless of their relationship.
    
    Args:
        definition: String containing the KEGG module definition
    
    Returns:
        List of KO IDs found in the definition
    """
    if not definition:
        return []
    
    # Use regex to extract all KO IDs (K followed by 5 digits)
    ko_pattern = r'K\d{5}'
    ko_ids = re.findall(ko_pattern, definition)
    
    # Return unique KO IDs
    return list(set(ko_ids))

def extract_custom_module_genes(definition):
    """Extract KO IDs from a custom module definition.
    
    Custom module definitions use a complex notation with +, -, and parentheses.
    This function extracts all KO IDs regardless of their relationship.
    
    Args:
        definition: String containing the custom module definition
    
    Returns:
        List of KO IDs found in the definition
    """
    if not definition:
        return []
    
    # Use regex to extract all KO IDs (K followed by 5 digits)
    ko_pattern = r'K\d{5}'
    ko_ids = re.findall(ko_pattern, definition)
    
    # Return unique KO IDs
    return list(set(ko_ids))

def load_module_data():
    """Load and parse module data from TSV file."""
    modules = {}
    custom_modules = {}
    
    # Add debug logging
    logging.info(f"Loading module data from {module_file}")
    
    try:
        # Read first few lines to debug
        with open(module_file, 'r') as f:
            first_lines = '\n'.join([f.readline() for _ in range(5)])
            logging.info(f"First lines of {module_file}:\n{first_lines}")
        
        # Now read the actual file
        with open(module_file, 'r') as f:
            # Skip the header row
            header = f.readline().strip()
            
            if 'Enceladus' in header and 'Europa' in header:
                logging.info("Found Enceladus and Europa columns in the header")
            else:
                logging.warning(f"Did not find Enceladus/Europa columns in header: {header}")
            
            for line in f:
                parts = line.strip().split('\t')
                
                # New format with Enceladus and Europa scores
                if len(parts) >= 4:  # At least 4 fields: Enceladus, Europa, Module ID, Description
                    try:
                        enceladus_score = int(parts[0]) if parts[0].strip() else 0
                        europa_score = int(parts[1]) if parts[1].strip() else 0
                        module_id = parts[2].strip()
                        description = parts[3].strip() if len(parts) > 3 else ""
                        
                        # Debug the extraction
                        if module_id:
                            logging.info(f"Parsed module: {module_id}, Enceladus: {enceladus_score}, Europa: {europa_score}")
                    except (ValueError, IndexError) as e:
                        logging.error(f"Error parsing line: {line}, Error: {e}")
                        continue
                else:
                    logging.warning(f"Skipping line with not enough fields: {line}")
                    continue
                
                # Extract genes based on the module description
                if module_id.startswith('M'):
                    # KEGG module
                    genes = extract_kegg_module_genes(description)
                    modules[module_id] = {
                        'genes': genes,
                        'description': description,
                        'name': description.split(' - ')[0] if ' - ' in description else description,
                        'enceladus_score': enceladus_score,
                        'europa_score': europa_score
                    }
                else:
                    # Custom module
                    genes = extract_custom_module_genes(description)
                    custom_modules[module_id] = {
                        'genes': genes,
                        'description': description,
                        'name': module_id.replace('_', ' '),
                        'enceladus_score': enceladus_score,
                        'europa_score': europa_score
                    }
    except Exception as e:
        logging.error(f"Error loading module data: {str(e)}")
        return {}, {}
    
    logging.info(f"Loaded {len(modules)} KEGG modules and {len(custom_modules)} custom modules")
    return modules, custom_modules

def calculate_module_completeness(hits_df, modules, custom_modules):
    """Calculate module completeness based on identified genes.
    
    For KEGG modules (starting with 'M'), calculate pairwise completeness similar to anvi_estimate_function.
    For custom modules, use adjusted frequency calculation: (found_genes/total_genes) * 1.5, capped at 100%.
    """
    logging.info("Calculating module completeness")
    
    # Process sequence results to extract gene IDs
    gene_ids, total_orfs, gene_ids_to_query_ids = process_sequence_results(hits_df)
    
    # Calculate gene-specific hits per million ORFs
    gene_hits_per_million = {}
    for gene_id, query_ids in gene_ids_to_query_ids.items():
        gene_hits_per_million[gene_id] = (len(query_ids) / total_orfs) * 1000000 if total_orfs > 0 else 0
    
    # Initialize results dictionary
    results = {}
    
    # Process standard KEGG modules (pairwise completeness)
    for module_id, module in modules.items():
        if not module.get('genes'):
            continue
            
        # Find genes present in the module
        found_genes = [gene for gene in module['genes'] if gene in gene_ids]
        total_genes = len(module['genes'])
        
        # Calculate completeness percentage using pairwise method for KEGG modules
        completeness = (len(found_genes) / total_genes) if total_genes > 0 else 0
        
        # Calculate normalized count (hits per million ORFs)
        count_per_million = 0
        if total_orfs > 0 and found_genes:
            # Count all query IDs (ORFs) that matched to any gene in this module
            module_orfs = set()
            for gene in found_genes:
                if gene in gene_ids_to_query_ids:
                    module_orfs.update(gene_ids_to_query_ids[gene])
            count_per_million = (len(module_orfs) / total_orfs) * 1000000
        
        # Create gene-specific data
        gene_data = {}
        for gene in module['genes']:
            if gene in gene_ids and gene in gene_hits_per_million:
                gene_data[gene] = {
                    'gene_id': gene,
                    'hits_per_million': gene_hits_per_million[gene],
                    'query_count': len(gene_ids_to_query_ids.get(gene, []))
                }
        
        # Store results
        results[module_id] = {
            'name': module.get('name', module_id),
            'description': module.get('description', ''),
            'total_genes': total_genes,
            'found_genes': found_genes,
            'completeness': completeness,
            'enceladus_score': module.get('enceladus_score', 0),
            'europa_score': module.get('europa_score', 0),
            'count_per_million': count_per_million,
            'genes': gene_data
        }
    
    # Process custom modules (adjusted frequency)
    for module_id, module in custom_modules.items():
        if not module.get('genes'):
            continue
            
        # Find genes present in the module
        found_genes = [gene for gene in module['genes'] if gene in gene_ids]
        total_genes = len(module['genes'])
        
        # Calculate completeness percentage using adjusted frequency for custom modules
        # Multiply by 1.5 and cap at 100%
        basic_completeness = (len(found_genes) / total_genes) if total_genes > 0 else 0
        adjusted_completeness = min(basic_completeness * 1.5, 1.0)
        
        # Calculate normalized count (hits per million ORFs)
        count_per_million = 0
        if total_orfs > 0 and found_genes:
            # Count all query IDs (ORFs) that matched to any gene in this module
            module_orfs = set()
            for gene in found_genes:
                if gene in gene_ids_to_query_ids:
                    module_orfs.update(gene_ids_to_query_ids[gene])
            count_per_million = (len(module_orfs) / total_orfs) * 1000000
        
        # Create gene-specific data
        gene_data = {}
        for gene in module['genes']:
            if gene in gene_ids and gene in gene_hits_per_million:
                gene_data[gene] = {
                    'gene_id': gene,
                    'hits_per_million': gene_hits_per_million[gene],
                    'query_count': len(gene_ids_to_query_ids.get(gene, []))
                }
        
        # Store results
        results[module_id] = {
            'name': module.get('name', module_id),
            'description': module.get('description', ''),
            'total_genes': total_genes,
            'found_genes': found_genes,
            'completeness': adjusted_completeness,
            'enceladus_score': module.get('enceladus_score', 0),
            'europa_score': module.get('europa_score', 0),
            'count_per_million': count_per_million,
            'genes': gene_data
        }
    
    # Log summary of results
    complete_modules = sum(1 for module in results.values() if module['completeness'] >= 1.0)
    partial_modules = sum(1 for module in results.values() if 0 < module['completeness'] < 1.0)
    logging.info(f"Found {complete_modules} complete and {partial_modules} partial modules out of {len(results)} total")
    
    return results

def process_sequences(input_file, threads=4):
    """Process input sequences."""
    global progress
    try:
        progress.update({
            'status': 'running',
            'stage': 'prodigal',
            'stage_text': 'Running Prodigal gene prediction...',
            'start_time': time.time(),
            'current_step_start_time': time.time(),
            'error': None,
            'results': None,
            'done': False
        })

        # Run Prodigal
        protein_seqs, gene_count = run_prodigal(input_file)
        progress['gene_count'] = gene_count
        progress['stage_progress'] = 33

        if not protein_seqs:
            raise Exception("No genes found by Prodigal")

        # Run DIAMOND search
        progress.update({
            'stage': 'diamond',
            'stage_text': 'Running DIAMOND search...',
            'current_step_start_time': time.time()
        })

        hits_df = batch_diamond_search(protein_seqs, threads)
        progress['stage_progress'] = 66

        # Process results
        progress.update({
            'stage': 'processing_results',
            'stage_text': 'Processing results...',
            'current_step_start_time': time.time()
        })

        # Load modules data
        modules, custom_modules = load_module_data()

        # Calculate completeness
        results = calculate_module_completeness(hits_df, modules, custom_modules)
        
        # Get the total number of unique gene identifiers and genes with functions
        gene_ids, total_orfs, gene_ids_to_query_ids = process_sequence_results(hits_df)
        
        # Store gene_ids_to_query_ids for CSV generation
        global gene_module_mapping
        gene_module_mapping = {}
        
        # Map genes to modules
        for module_id, module in results.items():
            if module.get('found_genes'):
                for gene in module.get('found_genes'):
                    if gene not in gene_module_mapping:
                        gene_module_mapping[gene] = []
                    gene_module_mapping[gene].append(module_id)
        
        # Store results for visualization - properly using app context
        with app.app_context():
            g.analysis_results = results
            g.diamond_hits = hits_df
        
        # Always keep a local reference to the latest results for quick access
        app.latest_analysis_results = results
        app.latest_diamond_hits = hits_df
        
        # Save results to a file for persistence
        results_path = os.path.join(app.config['UPLOAD_FOLDER'], 'module_results.json')
        try:
            with open(results_path, 'w') as f:
                json.dump(results, f)
        except Exception as e:
            logging.error(f"Error saving results to file: {str(e)}")
        
        # Update progress with results
        progress.update({
            'status': 'complete',
            'stage': 'complete',
            'stage_text': 'Analysis complete!',
            'stage_progress': 100,
            'results': results,
            'matched_genes_count': len(gene_ids),
            'genes_with_functions_count': total_orfs,
            'gene_ids_to_query_ids': gene_ids_to_query_ids,
            'visualization_available': True,
            'done': True
        })

        return results

    except Exception as e:
        logging.error(f"Error in process_sequences: {str(e)}")
        traceback.print_exc()
        progress.update({
            'status': 'error',
            'error': str(e),
            'stage_text': f'Error: {str(e)}',
            'done': True
        })
        return None

def allowed_file(filename):
    """Check if the file is an allowed type."""
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in {'fasta', 'fa', 'fna'}

@app.route('/')
def index():
    """Render the main page."""
    return render_template('index.html')

@app.route('/progress')
def get_progress():
    """Get current progress."""
    global progress
    
    # Check if we have real results available (from current analysis or previous run)
    visualization_available = False
    
    # First check g object
    if hasattr(g, 'analysis_results') and g.analysis_results:
        visualization_available = True
        logging.info("Visualization available from g.analysis_results")
    # Then check if progress status is complete
    elif progress.get('status') == 'complete':
        visualization_available = True
        logging.info("Visualization available from complete status")
    # Then check app's backup reference
    elif hasattr(app, 'latest_analysis_results') and app.latest_analysis_results:
        visualization_available = True
        logging.info("Visualization available from app backup")
    # Then check for exported heatmap data
    elif os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heatmap_data.json')):
        visualization_available = True
        logging.info("Visualization available from exported heatmap data")
    # Finally check for saved results
    else:
        results_path = os.path.join(app.config['UPLOAD_FOLDER'], 'module_results.json')
        if os.path.exists(results_path):
            try:
                with open(results_path, 'r') as f:
                    saved_results = json.load(f)
                    if saved_results and len(saved_results) > 0:
                        visualization_available = True
                        logging.info("Visualization available from saved results")
            except Exception as e:
                logging.error(f"Error checking saved results: {e}")
    
    # FORCE the visualization to be available if we have any results at all
    if progress.get('results') and len(progress.get('results', {})) > 0:
        visualization_available = True
        logging.info("Forcing visualization available due to results in progress object")
    
    # Set the flag in the progress object
    progress['visualization_available'] = visualization_available
    # Always set status and ensure it's 'complete' if we have results
    if progress.get('results') and len(progress.get('results', {})) > 0:
        progress['status'] = 'complete'
    else:
        progress['status'] = progress.get('status', 'unknown')  # Ensure status is always set
    
    logging.info(f"Progress update: visualization_available={visualization_available}, status={progress.get('status')}")
    
    return jsonify(progress)

@app.route('/upload', methods=['POST'])
def upload_file():
    """Handle file upload."""
    global progress
    
    try:
        if 'file' not in request.files:
            return jsonify({'success': False, 'message': 'No file part in the request.'})

        file = request.files['file']
        if file.filename == '':
            return jsonify({'success': False, 'message': 'No file selected for uploading.'})

        if file and allowed_file(file.filename):
            # Save the file
            filename = secure_filename(file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(file_path)
            
            # Reset progress
            progress.clear()
            progress.update({
                'status': 'running',
                'stage': 'prodigal',
                'stage_text': 'Starting analysis...',
                'stage_progress': 0,
                'gene_count': 0,
                'error': None,
                'results': None,
                'start_time': time.time(),
                'current_step_start_time': time.time(),
                'done': False
            })
            
            # Start processing in background
            thread = Thread(target=process_sequences, args=(file_path, 4))
            thread.daemon = True
            thread.start()
            
            return jsonify({'success': True, 'message': 'File uploaded and processing started.'})
        else:
            return jsonify({'success': False, 'message': 'Allowed file types are .fasta, .fa, .fna'})
    except Exception as e:
        error_msg = f"Error processing upload: {str(e)}"
        logging.error(error_msg)
        progress.update({
            'error': error_msg,
            'status': 'error',
            'stage_text': error_msg,
            'done': True
        })
        return jsonify({'success': False, 'message': error_msg})

@app.route('/test_with_subset')
def test_with_subset():
    """Process a subset of sequences for testing."""
    global progress
    
    try:
        # Reset progress
        progress.clear()
        progress.update({
            'status': 'running',
            'stage': 'prodigal',
            'stage_text': 'Starting test...',
            'stage_progress': 0,
            'gene_count': 0,
            'error': None,
            'results': None,
            'start_time': time.time(),
            'current_step_start_time': time.time(),
            'done': False
        })
        
        # Look for test file
        test_file = os.path.join(app.config['UPLOAD_FOLDER'], 'test_subset.fna')
        if not os.path.exists(test_file):
            raise Exception("Test file not found in uploads directory")
        
        # Start processing in background
        thread = Thread(target=process_sequences, args=(test_file, 4))
        thread.daemon = True
        thread.start()
        
        return jsonify({'success': True, 'message': 'Test started'})
        
    except Exception as e:
        error_msg = f"Error starting test: {str(e)}"
        logging.error(error_msg)
        progress.update({
            'error': error_msg,
            'status': 'error',
            'stage_text': error_msg,
            'done': True
        })
        return jsonify({'success': False, 'message': error_msg})

@app.route('/download_csv')
def download_csv():
    """Generate and download a CSV summary of results."""
    global progress, gene_module_mapping
    
    # Create a CSV in memory
    csv_data = StringIO()
    writer = csv.writer(csv_data)
    
    # Write header
    writer.writerow(['unique_functions', 'Associated_Modules', 'Genes', 'Percent_Identity', 'Alignment_Length', 'E-value', 'Bit_Score'])
    
    # Get DIAMOND hits dataframe if available
    diamond_hits = None
    if hasattr(g, 'diamond_hits') and g.diamond_hits is not None:
        diamond_hits = g.diamond_hits
        logging.info("Using diamond hits from g object for CSV export")
    elif hasattr(app, 'latest_diamond_hits') and app.latest_diamond_hits is not None:
        diamond_hits = app.latest_diamond_hits
        logging.info("Using diamond hits from app backup for CSV export")
    
    if diamond_hits is not None:
        logging.info(f"Found {len(diamond_hits)} DIAMOND hits for CSV export")
    else:
        logging.warning("No DIAMOND hits available for CSV export")
    
    # Write data
    if 'gene_ids_to_query_ids' in progress and progress['gene_ids_to_query_ids']:
        for gene_id, query_ids in progress['gene_ids_to_query_ids'].items():
            if query_ids:  # Only include hits > 0
                modules = ','.join(gene_module_mapping.get(gene_id, ['NA']))
                genes = ','.join(query_ids)
                
                # Default values for DIAMOND data
                pident = 'NA'
                length = 'NA'
                evalue = 'NA'
                bitscore = 'NA'
                
                # Try to get DIAMOND data for this gene
                if diamond_hits is not None:
                    # Find the best hit for this gene
                    gene_hits = diamond_hits[diamond_hits['sseqid'].str.contains(gene_id, na=False)]
                    if not gene_hits.empty:
                        best_hit = gene_hits.sort_values('bitscore', ascending=False).iloc[0]
                        pident = str(best_hit.get('pident', 'NA'))
                        length = str(best_hit.get('length', 'NA'))
                        evalue = str(best_hit.get('evalue', 'NA'))
                        bitscore = str(best_hit.get('bitscore', 'NA'))
                
                writer.writerow([gene_id, modules, genes, pident, length, evalue, bitscore])
    
    # Prepare response
    output = make_response(csv_data.getvalue())
    output.headers["Content-Disposition"] = "attachment; filename=Enceladupedia_summary.csv"
    output.headers["Content-type"] = "text/csv"
    
    return output

@app.route('/direct_heatmap')
def direct_heatmap():
    """Render the direct heatmap page with module data."""
    # Log that this endpoint was called
    logging.info("===== DIRECT HEATMAP ENDPOINT CALLED =====")
    
    # Try to get analysis results from various sources
    results = None
    diamond_hits = None
    
    # First try to get from g object
    if hasattr(g, 'analysis_results') and g.analysis_results is not None:
        results = g.analysis_results
        logging.info("Using analysis results from g object")
        if hasattr(g, 'diamond_hits') and g.diamond_hits is not None:
            diamond_hits = g.diamond_hits
    
    # Then try the app's backup reference
    elif hasattr(app, 'latest_analysis_results') and app.latest_analysis_results is not None:
        results = app.latest_analysis_results
        logging.info("Using analysis results from app backup")
        if hasattr(app, 'latest_diamond_hits') and app.latest_diamond_hits is not None:
            diamond_hits = app.latest_diamond_hits
    
    # Then try to load from saved file
    else:
        results_path = os.path.join(app.config['UPLOAD_FOLDER'], 'module_results.json')
        if os.path.exists(results_path):
            try:
                with open(results_path, 'r') as f:
                    results = json.load(f)
                    logging.info("Using analysis results from saved file")
                
                # Try to load diamond hits from saved file
                diamond_path = os.path.join(app.config['UPLOAD_FOLDER'], 'diamond_hits.csv')
                if os.path.exists(diamond_path):
                    try:
                        diamond_hits = pd.read_csv(diamond_path)
                        logging.info("Loaded diamond hits from saved file")
                    except Exception as e:
                        logging.error(f"Error loading diamond hits: {e}")
            except Exception as e:
                logging.error(f"Error loading saved results file: {e}")
    
    # As a last resort, try the exported data
    if not results:
        logging.info("No analysis results found, checking exported data")
        exported_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heatmap_data.json')
        if os.path.exists(exported_path):
            try:
                with open(exported_path, 'r') as f:
                    data = json.load(f)
                    if 'modules' in data:
                        results = data.get('modules', {})
                        logging.info(f"Using exported heatmap data with {len(results)} modules")
                    else:
                        # The file might contain direct module data
                        results = data
                        logging.info(f"Using direct module data with {len(results)} modules")
            except Exception as e:
                logging.error(f"Error loading exported heatmap data: {e}")
    
    # If no results are available, redirect to the main page with an error message
    if not results:
        logging.warning("No results available for heatmap visualization")
        return render_template('index.html', error="No analysis results available for visualization. Please run an analysis first.")
    
    # Load module descriptions from the module_descriptions.tsv file
    module_descriptions = {}
    try:
        module_desc_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'module_descriptions.tsv')
        if os.path.exists(module_desc_path):
            with open(module_desc_path, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            module_id = parts[0]
                            description = parts[1] if len(parts) > 1 else ''
                            module_descriptions[module_id] = description
            logging.info(f"Loaded {len(module_descriptions)} module descriptions")
    except Exception as e:
        logging.error(f"Error loading module descriptions: {e}")
    
    # Also load module descriptions from module_descriptions.txt (for all 18 modules)
    try:
        root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        module_desc_txt_path = os.path.join(root_dir, 'module_descriptions.txt')
        if os.path.exists(module_desc_txt_path):
            with open(module_desc_txt_path, 'r') as f:
                # Skip header line
                next(f, None)
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:  # Enceladus, Europa, Module ID, Description
                            module_id = parts[2]
                            description = parts[3] if len(parts) > 3 else ''
                            module_descriptions[module_id] = description
                            
                            # If module is not in results, add a placeholder
                            if module_id not in results:
                                enceladus_score = float(parts[0]) if parts[0].replace('.', '', 1).isdigit() else 0
                                europa_score = float(parts[1]) if parts[1].replace('.', '', 1).isdigit() else 0
                                results[module_id] = {
                                    'description': description,
                                    'completeness': 0.0,  # Default completeness
                                    'enceladus_score': enceladus_score,
                                    'europa_score': europa_score,
                                    'count_per_million': 0.0,  # Default count
                                    'found_genes': [],
                                    'total_genes': 0,
                                    'genes': {}
                                }
            logging.info(f"Loaded additional module descriptions from module_descriptions.txt")
    except Exception as e:
        logging.error(f"Error loading module descriptions from txt file: {e}")
    
    # Enhance results with module descriptions and ensure gene-level data is properly structured
    for module_id, module_data in results.items():
        if isinstance(module_data, dict):
            # Add description if missing
            if 'description' not in module_data or not module_data['description']:
                module_data['description'] = module_descriptions.get(module_id, '')
            
            # Ensure gene-level data is properly structured
            if 'genes' not in module_data:
                module_data['genes'] = {}
            
            # If we have found_genes but no gene-level data, create gene entries
            if 'found_genes' in module_data and module_data['found_genes'] and 'count_per_million' in module_data:
                found_genes = module_data['found_genes']
                total_count = module_data.get('count_per_million', 0)
                
                # Only create gene entries if the genes dict is empty
                if not module_data['genes'] and total_count > 0 and len(found_genes) > 0:
                    # Try to load diamond hits to get gene-specific counts
                    if diamond_hits is not None:
                        # Process diamond hits to get gene-specific counts
                        gene_ids_to_query_ids = {}
                        for _, row in diamond_hits.iterrows():
                            subject_id = row['sseqid'] if 'sseqid' in row else row.get('subject_id', '')
                            query_id = row['qseqid'] if 'qseqid' in row else row.get('query_id', '')
                            
                            # Use the same gene ID extraction logic as in process_sequence_results
                            gene_id = None
                            
                            # Case 1: Format KXXXXX_organism:gene
                            if '_' in subject_id and subject_id.startswith('K'):
                                ko_id = subject_id.split('_')[0]
                                if len(ko_id) >= 6:  # KO IDs like K00001
                                    gene_id = ko_id
                            
                            # Case 2: Format with | separators (old format)
                            elif '|' in subject_id:
                                parts = subject_id.split('|')
                                if parts and len(parts) > 1:
                                    potential_gene = parts[-1]
                                    if potential_gene.startswith('K') and len(potential_gene) >= 6:
                                        gene_id = potential_gene
                                    else:
                                        # This might be a GenBank ID
                                        gene_id = potential_gene
                            
                            # Case 3: Direct KO IDs without separators
                            elif subject_id.startswith('K') and len(subject_id) >= 6 and not '_' in subject_id:
                                gene_id = subject_id
                            
                            # Case 4: GenBank ID or other format
                            elif (re.match(r'^[A-Z0-9]+\.[0-9]+$', subject_id) or  # Accession.version
                                  re.match(r'^[A-Z][A-Z][A-Z][0-9][0-9]_[0-9]+$', subject_id) or  # Locus tag
                                  re.match(r'^[A-Z][0-9][A-Z][0-9][0-9]_[0-9]+$', subject_id) or  # Locus tag variation
                                  re.match(r'^[A-Z][0-9][A-Z][A-Z][0-9]_[0-9]+$', subject_id)):   # Another locus tag variation
                                gene_id = subject_id
                            
                            # Case 5: Any other ID - assume it could be a valid ID
                            else:
                                gene_id = subject_id
                            
                            # Store the gene ID if found and it's in our module's found_genes
                            if gene_id and gene_id in found_genes:
                                if gene_id not in gene_ids_to_query_ids:
                                    gene_ids_to_query_ids[gene_id] = []
                                if query_id not in gene_ids_to_query_ids[gene_id]:
                                    gene_ids_to_query_ids[gene_id].append(query_id)
                        
                        # Calculate hits per million for each gene
                        total_orfs = diamond_hits['qseqid'].nunique()
                        if total_orfs > 0:
                            for gene_id, query_ids in gene_ids_to_query_ids.items():
                                hits_per_million = (len(query_ids) / total_orfs) * 1000000
                                module_data['genes'][gene_id] = {
                                    'gene_id': gene_id,
                                    'hits_per_million': hits_per_million,
                                    'query_count': len(query_ids)
                                }
                    
                    # If we couldn't get gene-specific counts from diamond hits,
                    # distribute counts evenly among found genes
                    if not module_data['genes']:
                        count_per_gene = total_count / len(found_genes)
                        for gene in found_genes:
                            module_data['genes'][gene] = {
                                'gene_id': gene,
                                'hits_per_million': count_per_gene,
                                'query_count': 0  # Unknown query count
                            }
    
    # Count modules with ≥50% completeness for the data summary
    high_completeness_count = 0
    for module_id, module_data in results.items():
        if isinstance(module_data, dict) and module_data.get('completeness', 0) >= 0.5:  # 50% completeness threshold
            high_completeness_count += 1
    
    logging.info(f"Rendering heatmap with {len(results)} modules, {high_completeness_count} with ≥50% completeness")
    
    # Pass the data to the template
    return render_template('direct_heatmap.html', 
                           module_data=json.dumps(results))

@app.route('/export_module_data')
def export_module_data():
    """Export the module data to a JSON file for the standalone heatmap"""
    # Try to get analysis results from various sources
    results = None
    diamond_hits = None
    
    # First try to get from g object
    if hasattr(g, 'analysis_results') and g.analysis_results is not None:
        results = g.analysis_results
        if hasattr(g, 'diamond_hits') and g.diamond_hits is not None:
            diamond_hits = g.diamond_hits
    
    # Then try the app's backup reference
    elif hasattr(app, 'latest_analysis_results') and app.latest_analysis_results is not None:
        results = app.latest_analysis_results
        if hasattr(app, 'latest_diamond_hits') and app.latest_diamond_hits is not None:
            diamond_hits = app.latest_diamond_hits
    
    # Then try to load from saved file
    else:
        results_path = os.path.join(app.config['UPLOAD_FOLDER'], 'module_results.json')
        if os.path.exists(results_path):
            try:
                with open(results_path, 'r') as f:
                    results = json.load(f)
                
                # Try to load diamond hits from saved file
                diamond_path = os.path.join(app.config['UPLOAD_FOLDER'], 'diamond_hits.csv')
                if os.path.exists(diamond_path):
                    try:
                        diamond_hits = pd.read_csv(diamond_path)
                    except Exception as e:
                        logging.error(f"Error loading diamond hits: {e}")
            except Exception as e:
                logging.error(f"Error loading saved results file: {e}")
    
    # If we have real results, use them
    if results:
        # Process the results to ensure gene-specific data is included
        for module_id, module_data in results.items():
            if isinstance(module_data, dict):
                # Ensure gene-level data is properly structured
                if 'genes' not in module_data:
                    module_data['genes'] = {}
                
                # If we have found_genes but no gene-level data, create gene entries
                if 'found_genes' in module_data and module_data['found_genes'] and 'count_per_million' in module_data:
                    found_genes = module_data['found_genes']
                    total_count = module_data.get('count_per_million', 0)
                    
                    # Only create gene entries if the genes dict is empty
                    if not module_data['genes'] and total_count > 0 and len(found_genes) > 0:
                        # Try to get gene-specific counts from diamond hits
                        if diamond_hits is not None:
                            # Process diamond hits to get gene-specific counts
                            gene_ids_to_query_ids = {}
                            for _, row in diamond_hits.iterrows():
                                subject_id = row['sseqid'] if 'sseqid' in row else row.get('subject_id', '')
                                query_id = row['qseqid'] if 'qseqid' in row else row.get('query_id', '')
                                
                                # Use the same gene ID extraction logic as in process_sequence_results
                                gene_id = None
                                
                                # Case 1: Format KXXXXX_organism:gene
                                if '_' in subject_id and subject_id.startswith('K'):
                                    ko_id = subject_id.split('_')[0]
                                    if len(ko_id) >= 6:  # KO IDs like K00001
                                        gene_id = ko_id
                                
                                # Case 2: Format with | separators (old format)
                                elif '|' in subject_id:
                                    parts = subject_id.split('|')
                                    if parts and len(parts) > 1:
                                        potential_gene = parts[-1]
                                        if potential_gene.startswith('K') and len(potential_gene) >= 6:
                                            gene_id = potential_gene
                                        else:
                                            # This might be a GenBank ID
                                            gene_id = potential_gene
                                
                                # Case 3: Direct KO IDs without separators
                                elif subject_id.startswith('K') and len(subject_id) >= 6 and not '_' in subject_id:
                                    gene_id = subject_id
                                
                                # Case 4: GenBank ID or other format
                                elif (re.match(r'^[A-Z0-9]+\.[0-9]+$', subject_id) or  # Accession.version
                                      re.match(r'^[A-Z][A-Z][A-Z][0-9][0-9]_[0-9]+$', subject_id) or  # Locus tag
                                      re.match(r'^[A-Z][0-9][A-Z][0-9][0-9]_[0-9]+$', subject_id) or  # Locus tag variation
                                      re.match(r'^[A-Z][0-9][A-Z][A-Z][0-9]_[0-9]+$', subject_id)):   # Another locus tag variation
                                    gene_id = subject_id
                                
                                # Case 5: Any other ID - assume it could be a valid ID
                                else:
                                    gene_id = subject_id
                                
                                # Store the gene ID if found and it's in our module's found_genes
                                if gene_id and gene_id in found_genes:
                                    if gene_id not in gene_ids_to_query_ids:
                                        gene_ids_to_query_ids[gene_id] = []
                                    if query_id not in gene_ids_to_query_ids[gene_id]:
                                        gene_ids_to_query_ids[gene_id].append(query_id)
                            
                            # Calculate hits per million for each gene
                            total_orfs = diamond_hits['qseqid'].nunique()
                            if total_orfs > 0:
                                for gene_id, query_ids in gene_ids_to_query_ids.items():
                                    hits_per_million = (len(query_ids) / total_orfs) * 1000000
                                    module_data['genes'][gene_id] = {
                                        'gene_id': gene_id,
                                        'hits_per_million': hits_per_million,
                                        'query_count': len(query_ids)
                                    }
                        
                        # If we couldn't get gene-specific counts from diamond hits,
                        # distribute counts evenly among found genes
                        if not module_data['genes']:
                            count_per_gene = total_count / len(found_genes)
                            for gene in found_genes:
                                module_data['genes'][gene] = {
                                    'gene_id': gene,
                                    'hits_per_million': count_per_gene,
                                    'query_count': 0  # Unknown query count
                                }
        
        # Save to a JSON file
        output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heatmap_data.json')
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        return jsonify({
            'success': True,
            'message': f"Exported {len(results)} modules to heatmap_data.json"
        })
    else:
        # If no real results, load module data and create demo data
        modules, custom_modules = load_module_data()
        
        # Process modules to get completeness data
        results = {}
        
        # Process standard KEGG modules
        for module_id, module in modules.items():
            # Create gene-level data structure for the heatmap
            gene_data = {}
            module_genes = module.get('genes', [])
            for gene_id in module_genes:
                gene_data[gene_id] = {
                    'gene_id': gene_id,
                    'hits_per_million': random.uniform(10, 100),  # Random value for demonstration
                    'query_count': random.randint(1, 10)  # Random query count for demonstration
                }
                
            # Create a subset of genes as "found" for demonstration
            found_genes = random.sample(module_genes, min(random.randint(0, 3), len(module_genes)))
            
            results[module_id] = {
                'name': module.get('name', module_id),
                'description': module.get('description', ''),
                'total_genes': len(module_genes),
                'found_genes': found_genes,
                'completeness': random.uniform(0.1, 1.0), # Random completeness for visualization
                'enceladus_score': module.get('enceladus_score', 0),
                'europa_score': module.get('europa_score', 0),
                'count_per_million': random.uniform(20, 80),  # Random value for visualization
                'genes': gene_data  # Add gene-level data for the multi-section layout
            }
        
        # Process custom modules
        for module_id, module in custom_modules.items():
            # Create gene-level data structure for the heatmap
            gene_data = {}
            module_genes = module.get('genes', [])
            for gene_id in module_genes:
                gene_data[gene_id] = {
                    'gene_id': gene_id,
                    'hits_per_million': random.uniform(10, 100),  # Random value for demonstration
                    'query_count': random.randint(1, 10)  # Random query count for demonstration
                }
            
            # Create a subset of genes as "found" for demonstration
            found_genes = random.sample(module_genes, min(random.randint(0, 3), len(module_genes)))
            
            results[module_id] = {
                'name': module.get('name', module_id),
                'description': module.get('description', ''),
                'total_genes': len(module_genes),
                'found_genes': found_genes,
                'completeness': random.uniform(0.1, 1.0), # Random completeness for visualization
                'enceladus_score': module.get('enceladus_score', 0),
                'europa_score': module.get('europa_score', 0),
                'count_per_million': random.uniform(20, 80),  # Random value for visualization
                'genes': gene_data  # Add gene-level data for the multi-section layout
            }
        
        # Save to a JSON file
        output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heatmap_data.json')
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        return jsonify({
            'success': True,
            'message': f"Created demo data with {len(results)} modules in heatmap_data.json"
        })

@app.route('/export_heatmap_data', methods=['POST'])
def export_heatmap_data():
    """Export the current analysis results to a JSON file for the heatmap visualization."""
    try:
        # Try to get analysis results from various sources
        results = None
        
        # First try to get from g object
        if hasattr(g, 'analysis_results') and g.analysis_results is not None:
            results = g.analysis_results
            logging.info("Exporting results from g object")
        
        # Then try the app's backup reference
        elif hasattr(app, 'latest_analysis_results') and app.latest_analysis_results is not None:
            results = app.latest_analysis_results
            logging.info("Exporting results from app backup")
        
        # Then try to load from saved file
        else:
            results_path = os.path.join(app.config['UPLOAD_FOLDER'], 'module_results.json')
            if os.path.exists(results_path):
                try:
                    with open(results_path, 'r') as f:
                        results = json.load(f)
                        logging.info("Exporting results from saved file")
                except Exception as e:
                    logging.error(f"Error loading saved results file for export: {e}")
        
        if not results:
            return jsonify({"success": False, "message": "No analysis results available to export"})
        
        # Save the data to a file
        try:
            export_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heatmap_data.json')
            with open(export_path, 'w') as f:
                # Export the data directly without additional formatting
                # The standalone HTML will handle the data format
                json.dump(results, f, indent=2)
            
            logging.info(f"Heatmap data exported to {export_path}")
            return jsonify({'success': True, 'message': 'Heatmap data exported successfully'})
        except Exception as e:
            logging.error(f"Error exporting heatmap data: {e}")
            return jsonify({'error': f'Failed to export heatmap data: {str(e)}'}), 500
    except Exception as e:
        logging.error(f"Error exporting heatmap data: {e}")
        return jsonify({"success": False, "message": f"Error: {str(e)}"})

if __name__ == '__main__':
    app.run(debug=True, port=8888)