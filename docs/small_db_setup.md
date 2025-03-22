# Small Database Setup Guide

If you can't download the full UniProt-KEGG database (~20GB), you can create a smaller database focused on methanogenesis and related pathways.

## Option 1: Use Pre-filtered KEGG IDs

1. Create a list of relevant KEGG IDs:
```bash
# In your project root
mkdir -p databases/kegg
cd databases/kegg

# Create a file with your KEGG IDs
cat > kegg_ids.txt << EOL
K00399  # methyl-coenzyme M reductase alpha subunit
K00401  # methyl-coenzyme M reductase beta subunit
K00402  # methyl-coenzyme M reductase gamma subunit
# Add more KEGG IDs here
EOL
```

2. Download only relevant sequences:
```python
import requests
import time

def download_kegg_sequences(kegg_ids, output_file):
    """Download sequences for specific KEGG IDs."""
    with open(output_file, 'w') as f:
        for kid in kegg_ids:
            # Get protein sequences
            response = requests.get(f"https://rest.kegg.jp/get/{kid}/aaseq")
            if response.ok:
                f.write(f">{kid}\n{response.text.strip()}\n")
            time.sleep(1)  # Be nice to KEGG API

# Usage
with open('kegg_ids.txt') as f:
    kegg_ids = [line.split('#')[0].strip() for line in f if line.strip() and not line.startswith('#')]

download_kegg_sequences(kegg_ids, 'small_kegg.fasta')
```

3. Create DIAMOND database:
```bash
diamond makedb --in small_kegg.fasta --db small_kegg
```

4. Update app configuration:
Edit `app.py` to use the smaller database:
```python
app.config['DIAMOND_DB'] = 'databases/kegg/small_kegg'
```

## Option 2: Use Pre-computed Database

1. Download our pre-computed database:
```bash
# Replace URL with actual download link
curl -O https://your-server.com/small_kegg_db.dmnd
```

2. Place it in the correct directory:
```bash
mv small_kegg_db.dmnd databases/kegg/
```

3. Update app configuration:
```python
app.config['DIAMOND_DB'] = 'databases/kegg/small_kegg_db'
```

## Performance Comparison

| Database      | Size  | Setup Time | Search Speed | Coverage |
|--------------|-------|------------|--------------|----------|
| Full UniProt | ~20GB | 1-2 hours  | Fast        | Complete |
| Small DB     | ~500MB| 5 minutes  | Very Fast   | Limited  |

The small database is perfect for:
- Testing and development
- Quick analysis of specific pathways
- Systems with limited storage

Note: The small database only contains sequences related to methanogenesis and related pathways. If you need to search for other genes, use the full database.
