import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq # Added for reverse complement
import primer3
import os
import csv
import subprocess
import shutil
import multiprocessing
from functools import partial
from tqdm import tqdm

# --- Hardcoded Adapter Tails ---
FWD_TAIL = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
REV_TAIL = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# --- Core Helper Functions ---

def read_lines_from_file(file_path):
    """Reads a list of items from a file, one per line, skipping empty lines."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

# --- Mode 1: Design Pipeline Functions ---

def create_blast_db_if_needed(genome_fasta_path, blast_db_prefix):
    """
    Checks if a BLAST database exists. If not, it creates one from the provided
    genome FASTA, using Biopython to robustly clean headers and format.
    """
    db_files_exist = os.path.exists(f"{blast_db_prefix}.nin") or os.path.exists(f"{blast_db_prefix}.nhr")
    
    if db_files_exist:
        print(f"BLAST database '{blast_db_prefix}' already exists. Skipping creation.")
        return

    print(f"BLAST database not found. Creating a new one from '{genome_fasta_path}'...")
    
    if not shutil.which("makeblastdb"):
        raise RuntimeError("`makeblastdb` command not found. Please ensure NCBI BLAST+ is installed and in your system's PATH.")

    cleaned_fasta_path = f"{os.path.splitext(genome_fasta_path)[0]}.cleaned.fna"
    print(f"Using Biopython to create a clean FASTA file: '{cleaned_fasta_path}'")
    
    cleaned_records = []
    try:
        for record in SeqIO.parse(genome_fasta_path, "fasta"):
            record.id = record.id.split()[0]
            record.description = '' # Clear the description
            cleaned_records.append(record)
    except Exception as e:
         raise ValueError(f"Biopython could not parse '{genome_fasta_path}'. It may not be a valid FASTA file. Original error: {e}")

    if not cleaned_records:
        raise ValueError(f"Input FASTA file '{genome_fasta_path}' is empty or contains no valid records. Please check the file format.")

    SeqIO.write(cleaned_records, cleaned_fasta_path, "fasta")
    
    print(f"Building BLAST database '{blast_db_prefix}'...")
    command = ['makeblastdb', '-in', cleaned_fasta_path, '-dbtype', 'nucl', '-out', blast_db_prefix]
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
        print("BLAST database created successfully.")
    except subprocess.CalledProcessError as e:
        print("--- MAKEBLASTDB FAILED ---")
        print(f"Error: {e.stderr}")
        raise

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
        core_contig_id = contig_id.split('|')[-1] if '|' in contig_id else contig_id
        if core_contig_id not in genome_records:
             print(f"Warning: Contig '{contig_id}' (or '{core_contig_id}') for gene '{target_id}' not found in FASTA file. Skipping.")
             return None, None
        target_info['contig'] = core_contig_id
    
    contig_seq = genome_records[target_info['contig']].seq
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_NUM_RETURN': 20
    }
    return primer3.design_primers(seq_args, global_args)

def run_blast_specificity_check(primer_seq, blast_db_path):
    """Runs blastn for a single primer sequence and returns the number of exact matches."""
    command = [
        'blastn', '-query', '-', '-db', blast_db_path, '-task', 'blastn-short',
        '-outfmt', '6', '-perc_identity', '100', '-qcov_hsp_perc', '100'
    ]
    process = subprocess.run(command, input=f">primer\n{primer_seq}", capture_output=True, text=True, check=True)
    return len(process.stdout.strip().split('\n')) if process.stdout.strip() else 0

def write_design_output_files(all_csv_rows, all_bed_rows, output_prefix):
    """Writes all collected primer data from the design mode to consolidated files."""
    if not all_csv_rows:
        print("\nNo specific primers were successfully designed for any target.")
        return
    csv_file, bed_file = f"{output_prefix}.csv", f"{output_prefix}.bed"
    csv_headers = [
        'target_id', 'pair_rank', 
        'fwd_primer_tailed', 'rev_primer_tailed',
        'fwd_primer_seq', 'rev_primer_seq',
        'fwd_primer_tm', 'rev_primer_tm', 
        'amplicon_size', 'specificity_hits'
    ]
    with open(csv_file, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=csv_headers)
        writer.writeheader()
        writer.writerows(all_csv_rows)
    with open(bed_file, 'w') as bedf:
        bedf.writelines(all_bed_rows)
    print(f"\nResults for {len(all_bed_rows)} specific targets saved to '{csv_file}' and '{bed_file}'")

def run_design_mode(args):
    """Runs the full primer design and specificity check pipeline."""
    all_csv_rows, all_bed_rows = [], []
    
    create_blast_db_if_needed(args.genome, args.blast_db)
    
    cleaned_genome_path = f"{os.path.splitext(args.genome)[0]}.cleaned.fna"
    print(f"Loading CLEANED genome from '{cleaned_genome_path}'...")
    genome_records = SeqIO.to_dict(SeqIO.parse(cleaned_genome_path, "fasta"))
    
    print(f"Parsing GFF file '{args.gff}'...")
    gene_coords = parse_gff(args.gff)
    
    print(f"Reading target IDs from '{args.target_file}'...")
    target_ids = read_lines_from_file(args.target_file)
    
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

            if fwd_hits == 1 and rev_hits == 1:
                print(f"  -> SUCCESS: Found specific pair (Rank {i}).")
                
                fwd_rc = str(Seq(fwd_seq).reverse_complement())
                rev_rc = str(Seq(rev_seq).reverse_complement())
                
                fwd_tailed = fwd_rc + FWD_TAIL
                rev_tailed = rev_rc + REV_TAIL

                all_csv_rows.append({
                    'target_id': target_id, 'pair_rank': i, 
                    'fwd_primer_tailed': fwd_tailed,
                    'rev_primer_tailed': rev_tailed,
                    'fwd_primer_seq': fwd_seq, 
                    'rev_primer_seq': rev_seq,
                    'fwd_primer_tm': f"{primer_results[f'PRIMER_LEFT_{i}_TM']:.2f}",
                    'rev_primer_tm': f"{primer_results[f'PRIMER_RIGHT_{i}_TM']:.2f}",
                    'amplicon_size': primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                    'specificity_hits': f"F:{fwd_hits}, R:{rev_hits}"
                })
                
                product_size = primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                fwd_start_in_gene = primer_results[f'PRIMER_LEFT_{i}'][0]
                amp_start = target_info['start'] - 1 + fwd_start_in_gene
                amp_end = amp_start + product_size
                bed_row = f"{target_info['contig']}\t{amp_start}\t{amp_end}\t{target_id}_amplicon_{i}\t0\t+\n"
                all_bed_rows.append(bed_row)
                
                specific_pair_found = True
                break 
        
        if not specific_pair_found:
             print(f"  -> Could not find a specific primer pair for {target_id} among the candidates.")

    write_design_output_files(all_csv_rows, all_bed_rows, args.output_prefix)

# --- Mode 2: Tail-Only Pipeline Functions ---

def write_tail_only_csv(all_csv_rows, output_prefix):
    """Writes the results of the tail-only mode to a CSV file."""
    if not all_csv_rows:
        print("\nNo primers were processed.")
        return
    
    csv_file = f"{output_prefix}.csv"
    csv_headers = [
        'pair_id', 'fwd_primer_tailed', 'rev_primer_tailed',
        'fwd_primer_seq', 'rev_primer_seq'
    ]
    with open(csv_file, 'w', newline='') as csvf:
        writer = csv.DictWriter(csvf, fieldnames=csv_headers)
        writer.writeheader()
        writer.writerows(all_csv_rows)
    print(f"\nSuccessfully tailed {len(all_csv_rows)} primer pairs.")
    print(f"Results saved to '{csv_file}'")

def run_tail_only_mode(args):
    """Runs the logic to simply tail existing primer files."""
    print(f"Reading forward primers from: {args.tail_fwd_file}")
    print(f"Reading reverse primers from: {args.tail_rev_file}")
    
    fwd_primers = read_lines_from_file(args.tail_fwd_file)
    rev_primers = read_lines_from_file(args.tail_rev_file)
    
    if len(fwd_primers) != len(rev_primers):
        print("\nWarning: The number of forward and reverse primers does not match.")
        print(f"  Forward primers found: {len(fwd_primers)}")
        print(f"  Reverse primers found: {len(rev_primers)}")
        print("Processing the minimum number of pairs.")
    
    all_csv_rows = []
    num_pairs = min(len(fwd_primers), len(rev_primers))
    
    for i in range(num_pairs):
        fwd_seq = fwd_primers[i]
        rev_seq = rev_primers[i]
        
        fwd_rc = str(Seq(fwd_seq).reverse_complement())
        rev_rc = str(Seq(rev_seq).reverse_complement())
        
        fwd_tailed = fwd_rc + FWD_TAIL
        rev_tailed = rev_rc + REV_TAIL
        
        all_csv_rows.append({
            'pair_id': f"pair_{i+1}",
            'fwd_primer_tailed': fwd_tailed,
            'rev_primer_tailed': rev_tailed,
            'fwd_primer_seq': fwd_seq,
            'rev_primer_seq': rev_seq
        })
        
    write_tail_only_csv(all_csv_rows, args.output_prefix)

# --- Main Execution "Router" ---

def main():
    parser = argparse.ArgumentParser(description="Design or tail multiplex PCR primers.")
    
    # --- Design Mode Arguments ---
    parser.add_argument('--genome', help="Path to the reference genome in FASTA format.")
    parser.add_argument('--gff', help="Path to a GFF file for gene coordinate lookups.")
    parser.add_argument('--target-file', help="Path to a text file with one target gene ID per line.")
    parser.add_argument('--blast-db', help="Prefix for the local BLAST database.")
    
    # --- Tail-Only Mode Arguments ---
    parser.add_argument('--tail-fwd-file', help="Path to a file of forward primers (one per line).")
    parser.add_argument('--tail-rev-file', help="Path to a file of reverse primers (one per line).")
    
    # --- Common Argument ---
    parser.add_argument('--output-prefix', default='final_primers', help="Prefix for output files.")
    
    args = parser.parse_args()

    try:
        # Check which mode to run
        if args.tail_fwd_file and args.tail_rev_file:
            print("--- Running in Tail-Only Mode ---")
            run_tail_only_mode(args)
            
        elif args.genome and args.gff and args.target_file and args.blast_db:
            print("--- Running in Full Design Mode ---")
            run_design_mode(args)
            
        else:
            print("Error: You must provide arguments for either 'Design Mode' or 'Tail-Only Mode'.")
            print("\nFor Design Mode, specify: --genome, --gff, --target-file, and --blast-db")
            print("For Tail-Only Mode, specify: --tail-fwd-file and --tail-rev-file")

    except (FileNotFoundError, subprocess.CalledProcessError, RuntimeError, ValueError) as e:
        print(f"\nAn error occurred: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    main()

