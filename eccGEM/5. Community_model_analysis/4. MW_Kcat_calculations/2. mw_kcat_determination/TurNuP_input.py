import cobra
import pandas as pd
import re
from Bio import SeqIO
from collections import defaultdict

def load_fasta_sequences(fasta_files):
    """Load sequences from FASTA files"""
    sequences = {}
    for fasta_file in fasta_files:
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences[record.id] = str(record.seq)
        except FileNotFoundError:
            print(f"Warning: {fasta_file} not found")
    return sequences

def extract_genes_from_gpr(gpr_string):
    """Extract all individual genes from GPR string"""
    if not gpr_string or gpr_string.strip() == "":
        return []
    
    # Remove parentheses and boolean operators, extract gene names
    # This regex finds sequences of word characters (letters, digits, underscores)
    gene_pattern = r'\b[A-Za-z0-9_]+\b'
    potential_genes = re.findall(gene_pattern, gpr_string)
    
    # Filter out boolean operators
    boolean_ops = {'and', 'or', 'AND', 'OR', 'not', 'NOT'}
    genes = [gene for gene in potential_genes if gene not in boolean_ops]
    
    return list(set(genes))  # Remove duplicates

def get_kegg_id(metabolite):
    """Extract KEGG ID from metabolite annotations, removing suffixes"""
    
    def extract_clean_kegg_id(text):
        """Extract KEGG ID in format CXXXXX from text, removing suffixes"""
        if not text:
            return ""
        
        # Look for KEGG compound pattern (C followed by 5 digits)
        kegg_pattern = r'C\d{5}'
        match = re.search(kegg_pattern, str(text))
        if match:
            return match.group()
        return ""
    
    # Check metabolite annotations
    if hasattr(metabolite, 'annotation'):
        # Check for direct KEGG annotation
        if 'kegg.compound' in metabolite.annotation:
            kegg_ids = metabolite.annotation['kegg.compound']
            if isinstance(kegg_ids, list):
                for kegg_id in kegg_ids:
                    clean_id = extract_clean_kegg_id(kegg_id)
                    if clean_id:
                        return clean_id
            else:
                clean_id = extract_clean_kegg_id(kegg_ids)
                if clean_id:
                    return clean_id
        
        # Check for other KEGG annotation formats
        for key, value in metabolite.annotation.items():
            if 'kegg' in key.lower():
                if isinstance(value, list):
                    for v in value:
                        clean_id = extract_clean_kegg_id(v)
                        if clean_id:
                            return clean_id
                else:
                    clean_id = extract_clean_kegg_id(value)
                    if clean_id:
                        return clean_id
    
    # Check if metabolite ID itself contains a KEGG ID (like C00001_p)
    clean_id = extract_clean_kegg_id(metabolite.id)
    if clean_id:
        return clean_id
    
    # Check metabolite name
    if hasattr(metabolite, 'name'):
        clean_id = extract_clean_kegg_id(metabolite.name)
        if clean_id:
            return clean_id
    
    return ""

# Load the model
print("Loading model...")
model = cobra.io.read_sbml_model("irreversible_cochliobolus_maize_leaf_community.xml")

# Load FASTA sequences
print("Loading FASTA sequences...")
fasta_files = ["CocheC4.fasta", "maize.fasta"]
sequences = load_fasta_sequences(fasta_files)

# Prepare data for Excel
data = []
skipped_genes = set()

print("Processing reactions...")
for reaction in model.reactions:
    # Only process reactions with GPR rules
    if not reaction.gene_reaction_rule or reaction.gene_reaction_rule.strip() == "":
        continue
    
    print(f"Processing reaction: {reaction.id}")
    
    # Extract genes from GPR
    genes = extract_genes_from_gpr(reaction.gene_reaction_rule)
    
    if not genes:
        continue
    
    # Get substrates (reactants)
    substrates = [met for met in reaction.reactants]
    
    # Create rows for each gene × substrate combination
    for gene in genes:
        # Only process genes that are available in FASTA files
        if gene in sequences:
            enzyme_sequence = sequences[gene]
            
            for substrate in substrates:
                # Get KEGG ID for substrate
                kegg_id = get_kegg_id(substrate)
                
                # Add row to data
                data.append({
                    'Reactions': reaction.id,
                    'Enzyme': enzyme_sequence,
                    'Substrates': kegg_id,
                    'Products': ""  # Always empty as requested
                })
        else:
            # Track skipped genes for reporting
            skipped_genes.add(gene)

# Create DataFrame and save to Excel
print("Creating Excel file...")
df = pd.DataFrame(data)

# Save to Excel
with pd.ExcelWriter('input.xlsx', engine='openpyxl') as writer:
    df.to_excel(writer, sheet_name='Sheet1', index=False)

print(f"Excel file created successfully!")
print(f"Total rows: {len(df)}")
print(f"Total reactions processed: {len(set(df['Reactions']))}")
print(f"Total unique genes found in FASTA: {len(sequences)}")
print(f"Total genes skipped (not in FASTA): {len(skipped_genes)}")

if skipped_genes:
    print(f"Skipped genes: {', '.join(sorted(skipped_genes))}")

# Print some examples for verification
print("\nFirst few rows preview:")
print(df.head())