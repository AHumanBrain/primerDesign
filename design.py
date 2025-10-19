import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import primer3
import os

# --- Core Functions ---

def parse_gff(gff_file):
    """
    Parses a GFF file to extract gene coordinates.
    Returns a dictionary mapping gene names to their coordinate info.
    """
    gene_coords = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            gene_name = None
            for attr in attributes.split(';'):
                if 'Name=' in attr:
                    gene_name = attr.split('=')[1]
                    break
            
            if gene_name:
                gene_coords[gene_name] = {
                    'contig': contig,
                    'start': start,
                    'end': end,
                    'strand': strand
                }
    return gene_coords

def extract_target_sequence(genome_records, gene_coords, target_id):
    """
    Extracts the DNA sequence for a target gene from the genome.
    """
    if target_id not in gene_coords:
        raise ValueError(f"Error: Target ID '{target_id}' not found in GFF file.")
    
    target_info = gene_coords[target_id]
    contig_id = target_info['contig']
    
    if contig_id not in genome_records:
        raise ValueError(f"Error: Contig '{contig_id}' for gene '{target_id}' not found in FASTA file.")
        
    contig_seq = genome_records[contig_id].seq
    
    # GFF coordinates are 1-based, Biopython is 0-based
    start = target_info['start'] - 1
    end = target_info['end']
    
    target_seq = contig_seq[start:end]
    
    return str(target_seq)

def design_primers_for_sequence(sequence, target_id):
    """
    Uses the primer3-py library to design primers for a given sequence.
    """
    # Primer3 expects a dictionary of arguments.
    # We define a SEQUENCE_ID and the sequence itself.
    seq_args = {
        'SEQUENCE_ID': target_id,
        'SEQUENCE_TEMPLATE': sequence,
    }
    
    # Global arguments for primer design criteria.
    # These will eventually come from a config file.
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]],
    }
    
    # The primer3.designPrimers function returns a dictionary of results.
    primer3_results = primer3.designPrimers(seq_args, global_args)
    
    return primer3_results

def print_primer_results(results):
    """
    Prints the designed primer pairs in a user-friendly format.
    """
    num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)
    print(f"--- Primer Design Results ---")
    if num_returned > 0:
        print(f"Successfully designed {num_returned} primer pair(s).\n")
        
        # We'll just show the best one (pair 0) for now.
        for i in range(min(1, num_returned)):
            print(f"Best Pair (Pair {i}):")
            fwd_seq = results[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev_seq = results[f'PRIMER_RIGHT_{i}_SEQUENCE']
            fwd_tm = results[f'PRIMER_LEFT_{i}_TM']
            rev_tm = results[f'PRIMER_RIGHT_{i}_TM']
            product_size = results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
            
            print(f"  Forward Primer: {fwd_seq}")
            print(f"  Reverse Primer: {rev_seq}")
            print(f"  Forward Tm:     {fwd_tm:.2f}°C")
            print(f"  Reverse Tm:     {rev_tm:.2f}°C")
            print(f"  Amplicon Size:  {product_size} bp")
            
    else:
        print("Could not find any suitable primers for the given target and parameters.")
        print("Reason:", results.get('PRIMER_EXPLAIN'))


# --- Main Execution Block ---

def main():
    """Main function to parse arguments and run the primer design workflow."""
    parser = argparse.ArgumentParser(description="Design multiplex PCR primers from a reference genome.")
    parser.add_argument('--genome', required=True, help="Path to the reference genome in FASTA format.")
    parser.add_argument('--gff', required=True, help="Path to a GFF file for gene coordinate lookups.")
    parser.add_argument('--target-id', required=True, help="The gene name or ID to design primers for.")
    # Placeholders for future functionality
    # parser.add_argument('--blast-db', required=True, help="Path to the local BLAST database of the genome.")
    # parser.add_argument('--output-prefix', default='primers', help="Prefix for output files.")
    
    args = parser.parse_args()

    # --- Workflow ---
    try:
        # 1. Load genome into memory
        print(f"Loading genome from '{args.genome}'...")
        genome_records = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
        
        # 2. Parse GFF to get gene coordinates
        print(f"Parsing GFF file '{args.gff}'...")
        gene_coords = parse_gff(args.gff)
        
        # 3. Extract the target sequence
        print(f"Extracting sequence for target '{args.target_id}'...")
        target_sequence = extract_target_sequence(genome_records, gene_coords, args.target_id)
        
        # 4. Design primers using Primer3
        print("Designing primers with Primer3...")
        primer_results = design_primers_for_sequence(target_sequence, args.target_id)
        
        # 5. Print the results
        print_primer_results(primer_results)

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
    except ValueError as e:
        print(e)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
