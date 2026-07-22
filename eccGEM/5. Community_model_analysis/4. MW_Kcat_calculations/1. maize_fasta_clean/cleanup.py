import cobra
from Bio import SeqIO

def extract_genes_from_model(model_file):
    """Extract gene IDs from COBRA model using cobrapy"""
    try:
        # Load the model
        model = cobra.io.read_sbml_model(model_file)
        
        # Extract all gene IDs from the model
        gene_ids = {gene.id for gene in model.genes}
        
        return gene_ids
    
    except Exception as e:
        print(f"Error loading model: {e}")
        return set()

def filter_fasta_with_biopython(input_fasta, output_fasta, gene_ids):
    """Filter FASTA file using BioPython to keep only genes present in gene_ids set"""
    sequences_kept = 0
    sequences_total = 0
    
    with open(output_fasta, 'w') as outfile:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequences_total += 1
            
            # Check if this gene ID is in our model
            if record.id in gene_ids:
                # Write header and sequence on single line
                outfile.write(f">{record.id}\n")
                outfile.write(f"{str(record.seq)}\n")
                sequences_kept += 1
    
    return sequences_kept, sequences_total

def filter_fasta_simple(input_fasta, output_fasta, gene_ids):
    """Filter FASTA file without BioPython dependency - single line sequences"""
    sequences_kept = 0
    sequences_total = 0
    
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        write_sequence = False
        current_header = ""
        current_sequence = ""
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                # Write previous sequence if it was selected
                if write_sequence and current_sequence:
                    outfile.write(f"{current_header}\n")
                    outfile.write(f"{current_sequence}\n")
                    sequences_kept += 1
                
                # Reset for new sequence
                sequences_total += 1
                gene_id = line[1:]  # Remove '>' character
                current_header = line
                current_sequence = ""
                
                # Check if this gene should be kept
                write_sequence = gene_id in gene_ids
                
            else:
                # This is sequence data - concatenate to single line
                if write_sequence:
                    current_sequence += line
        
        # Write the last sequence if it was selected
        if write_sequence and current_sequence:
            outfile.write(f"{current_header}\n")
            outfile.write(f"{current_sequence}\n")
            sequences_kept += 1
    
    return sequences_kept, sequences_total

def main():
    # File paths
    model_file = "izma6518_L_exchange_added.xml"
    input_fasta = "maize_draft.fasta"
    output_fasta = "maize.fasta"
    
    try:
        # Extract gene IDs from COBRA model
        print("Loading metabolic model and extracting genes...")
        gene_ids = extract_genes_from_model(model_file)
        print(f"Found {len(gene_ids)} genes in the metabolic model")
        
        if len(gene_ids) == 0:
            print("Warning: No genes found in the model")
            return
        
        # Print first few gene IDs for verification
        sample_genes = sorted(list(gene_ids))[:10]
        print("Sample gene IDs from model:", sample_genes)
        
        # Try to use BioPython for cleaner FASTA handling, fall back to simple method
        try:
            print("Filtering FASTA file using BioPython...")
            kept, total = filter_fasta_with_biopython(input_fasta, output_fasta, gene_ids)
        except ImportError:
            print("BioPython not available, using simple FASTA parser...")
            kept, total = filter_fasta_simple(input_fasta, output_fasta, gene_ids)
        
        print(f"\nFiltering complete!")
        print(f"Total sequences in input: {total}")
        print(f"Sequences kept: {kept}")
        print(f"Sequences removed: {total - kept}")
        print(f"Percentage kept: {(kept/total)*100:.1f}%")
        print(f"Filtered FASTA written to {output_fasta}")
        print("Each sequence is written on a single line")
        
        # Show which genes from the model were found in the FASTA
        if kept > 0:
            print(f"\nSuccessfully matched {kept} genes between model and FASTA file")
        else:
            print("\nWarning: No genes from the model were found in the FASTA file")
            print("This might indicate a mismatch in gene ID formatting")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        print("Please make sure both files are in the current directory")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()