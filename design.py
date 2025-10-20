import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import primer3
import os
import csv
import subprocess # New import for running external commands
import tempfile   # New import for creating temporary files

# --- Core Functions ---

def read_target_ids(file_path):
    """Reads a list of target IDs from a file, one per line."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def parse_gff(gff_file):
    """Parses a GFF file to extract gene coordinates."""
    gene_coords = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene': continue
            contig, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attributes = parts[8]
            gene_name = next((attr.split('=')[1] for attr in attributes.split(';') if 'Name=' in attr), None)
            if gene_name:
                gene_coords[gene_name] = {'contig': contig, 'start': start, 'end': end, 'strand': strand}
    return gene_coords

def extract_target_sequence(genome_records, gene_coords, target_id):
    """Extracts the DNA sequence for a target gene from the genome."""
    if target_id not in gene_coords:
        print(f"Warning: Target ID '{target_id}' not found in GFF file. Skipping.")
        return None, None
    target_info = gene_coords[target_id]
    contig_id = target_info['contig']
    if contig_id not in genome_records:
        print(f"Warning: Contig '{contig_id}' for gene '{target_id}' not found in FASTA file. Skipping.")
        return None, None
    contig_seq = genome_records[contig_id].seq
    start, end = target_info['start'] - 1, target_info['end']
    target_seq = contig_seq[start:end]
    return str(target_seq), target_info

def design_primers_for_sequence(sequence, target_id):
    """Uses primer3-py to design primers for a given sequence."""
    seq_args = {'SEQUENCE_ID': target_id, 'SEQUENCE_TEMPLATE': sequence}
    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0, 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0, 'PRIMER_MAX_GC': 80.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]],
        'PRIMER_NUM_RETURN': 20 # Ask for more primers to increase chance of finding a specific one
    }
    return primer3.design_primers(seq_args, global_args)

def run_blast_specificity_check(primer_seq, blast_db_path):
    """
    Runs blastn for a single primer sequence against a local BLAST database.
    Returns the number of exact matches found.
    """
    # blastn is very particular about short sequences, so we use specific parameters.
    # -task blastn-short is optimized for sequences < 30nt.
    # -outfmt 6 gives a simple, tab-delimited output.
    # -perc_identity 100 and -qcov_hsp_perc 100 ensure we only find perfect, full-length matches.
    command = [
        'blastn',
        '-query', '-',  # Read query from stdin
        '-db', blast_db_path,
        '-task', 'blastn-short',
        '-outfmt', '6',
        '-perc_identity', '100',
        '-qcov_hsp_perc', '100'
    ]
    
    # We pass the primer sequence to BLAST via standard input (stdin)
    process = subprocess.run(
        command,
        input=f">primer\n{primer_seq}",
        capture_output=True,
        text=True,
        check=True
    )
    
    # The number of lines in the output is the number of perfect hits.
    return len(process.stdout.strip().split('\n')) if process.stdout.strip() else 0

def write_output_files(all_csv_rows, all_bed_rows, output_prefix):
    """Writes all collected primer data to consolidated files."""
    if not all_csv_rows:
        print("\nNo specific primers were successfully designed for any target.")
        return

    csv_file = f"{output_prefix}.csv"
    bed_file = f"{output_prefix}.bed"
    
    csv_headers = [
        'target_id', 'pair_rank', 'fwd_primer_seq', 'fwd_primer_tm', 
        'rev_primer_seq', 'rev_primer_tm', 'amplicon_size', 'specificity_hits'
    ]
    with open(csv_file, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=csv_headers)
        writer.writeheader()
        writer.writerows(all_csv_rows)

    with open(bed_file, 'w') as bedf:
        bedf.writelines(all_bed_rows)

    print(f"\nResults for {len(all_bed_rows)} specific targets saved to '{csv_file}' and '{bed_file}'")

# --- Main Execution Block ---

def main():
    """Main function to parse arguments and run the primer design workflow."""
    parser = argparse.ArgumentParser(description="Design specific primers for multiple targets.")
    parser.add_argument('--genome', required=True, help="Path to the reference genome in FASTA format.")
    parser.add_argument('--gff', required=True, help="Path to a GFF file for gene coordinate lookups.")
    parser.add_argument('--target-file', required=True, help="Path to a text file with one target gene ID per line.")
    parser.add_argument('--blast-db', required=True, help="Path/name of the local BLAST database (e.g., 'ecoli_db').")
    parser.add_argument('--output-prefix', default='specific_multiplex_primers', help="Prefix for output files.")
    
    args = parser.parse_args()

    all_csv_rows = []
    all_bed_rows = []

    try:
        print(f"Loading genome from '{args.genome}'...")
        genome_records = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
        
        print(f"Parsing GFF file '{args.gff}'...")
        gene_coords = parse_gff(args.gff)
        
        print(f"Reading target IDs from '{args.target_file}'...") # <-- FIXED TYPO HERE
        target_ids = read_target_ids(args.target_file)
        
        for target_id in target_ids:
            print(f"\n--- Processing target: {target_id} ---")
            target_sequence, target_info = extract_target_sequence(genome_records, gene_coords, target_id)
            if not target_sequence: continue
            
            print("Designing candidate primers with Primer3...")
            primer_results = design_primers_for_sequence(target_sequence, target_id)
            
            num_returned = primer_results.get('PRIMER_PAIR_NUM_RETURNED', 0)
            if num_returned == 0:
                print("Primer3 could not design any initial primers for this target.")
                continue

            print(f"Found {num_returned} candidate pairs. Checking specificity with BLAST...")
            
            specific_pair_found = False
            for i in range(num_returned):
                fwd_seq = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
                rev_seq = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']

                fwd_hits = run_blast_specificity_check(fwd_seq, args.blast_db)
                rev_hits = run_blast_specificity_check(rev_seq, args.blast_db)

                # A primer pair is specific if BOTH primers hit the genome exactly ONCE.
                if fwd_hits == 1 and rev_hits == 1:
                    print(f"  -> Found specific pair (Rank {i}): Fwd hits: {fwd_hits}, Rev hits: {rev_hits}. Accepting.")
                    
                    # This is the first specific pair we've found for this target, so we take it.
                    all_csv_rows.append({
                        'target_id': target_id, 'pair_rank': i,
                        'fwd_primer_seq': fwd_seq,
                        'fwd_primer_tm': f"{primer_results[f'PRIMER_LEFT_{i}_TM']:.2f}",
                        'rev_primer_seq': rev_seq,
                        'rev_primer_tm': f"{primer_results[f'PRIMER_RIGHT_{i}_TM']:.2f}",
                        'amplicon_size': primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                        'specificity_hits': f"F:{fwd_hits}, R:{rev_hits}"
                    })
                    
                    product_size = primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                    fwd_primer_start_in_gene = primer_results[f'PRIMER_LEFT_{i}'][0]
                    amplicon_genome_start = target_info['start'] - 1 + fwd_primer_start_in_gene
                    amplicon_genome_end = amplicon_genome_start + product_size
                    bed_name = f"{target_id}_amplicon_{i}"
                    bed_row = (f"{target_info['contig']}\t{amplicon_genome_start}\t"
                               f"{amplicon_genome_end}\t{bed_name}\t0\t+\n")
                    all_bed_rows.append(bed_row)
                    
                    specific_pair_found = True
                    break # Move to the next target gene
            
            if not specific_pair_found:
                print(f"  -> Could not find a specific primer pair for {target_id} among the candidates.")

        write_output_files(all_csv_rows, all_bed_rows, args.output_prefix)

    except FileNotFoundError as e:
        print(f"Error: Input file not found - {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()

