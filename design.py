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
import itertools
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
        # Return a non-None error message
        return None, None, f"Warning: Target ID '{target_id}' not found in GFF file. Skipping."
    
    target_info = gene_coords[target_id]
    contig_id = target_info['contig']
    
    if contig_id not in genome_records:
        core_contig_id = contig_id.split('|')[-1] if '|' in contig_id else contig_id
        if core_contig_id not in genome_records:
             # Return a non-None error message
             return None, None, f"Warning: Contig '{contig_id}' (or '{core_contig_id}') for gene '{target_id}' not found in FASTA file. Skipping."
        target_info['contig'] = core_contig_id
    
    contig_seq = genome_records[target_info['contig']].seq
    start, end = target_info['start'] - 1, target_info['end']
    target_seq = contig_seq[start:end]
    # Return a None error message on success
    return str(target_seq), target_info, None

def design_primers_for_sequence(sequence, target_id, strategy_settings):
    """Uses primer3-py to design primers for a given sequence."""
    seq_args = {'SEQUENCE_ID': target_id, 'SEQUENCE_TEMPLATE': sequence}
    
    # Start with default global settings
    global_args = {
        'PRIMER_OPT_SIZE': 20, 'PRIMER_MIN_SIZE': 18, 'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0, 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0, 'PRIMER_MAX_GC': 80.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], # Default
        'PRIMER_NUM_RETURN': 20
    }
    
    # Apply strategy-specific settings, overriding defaults
    global_args.update(strategy_settings)
    
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

def process_single_target(target_id, genome_records, gene_coords, blast_db):
    """
    This function contains all the work for designing and checking primers
    for a SINGLE target. It finds ALL specific pairs across multiple strategies.
    Returns: (target_id, list_of_primer_data_dicts, error_message_or_None)
    """
    # 1. Define Retry Strategies
    strategies = [
        {'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}, # Strategy 0 (Default)
        {'PRIMER_PRODUCT_SIZE_RANGE': [[250, 350]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}, # Strategy 1 (Longer)
        {'PRIMER_PRODUCT_SIZE_RANGE': [[100, 150]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}, # Strategy 2 (Shorter)
        {'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_MIN_TM': 55.0, 'PRIMER_MAX_TM': 65.0}, # Strategy 3 (Relaxed Tm)
    ]
    
    # 2. Extract Sequence
    target_sequence, target_info, error = extract_target_sequence(genome_records, gene_coords, target_id)
    if not target_sequence:
        return (target_id, [], error) # Return target_id, empty list, and error

    all_specific_pairs_for_target = []

    # 3. Loop through all strategies
    for strategy_index, strategy_settings in enumerate(strategies):
        primer_results = design_primers_for_sequence(target_sequence, target_id, strategy_settings)
        num_returned = primer_results.get('PRIMER_PAIR_NUM_RETURNED', 0)
        if num_returned == 0:
            continue # Try next strategy

        # 4. Check Specificity for all pairs in this batch
        for i in range(num_returned):
            fwd_seq = primer_results[f'PRIMER_LEFT_{i}_SEQUENCE']
            rev_seq = primer_results[f'PRIMER_RIGHT_{i}_SEQUENCE']
            
            try:
                fwd_hits = run_blast_specificity_check(fwd_seq, blast_db)
                rev_hits = run_blast_specificity_check(rev_seq, blast_db)
            except Exception as e:
                return (target_id, [], f"BLAST failed for {target_id}: {e}")

            if fwd_hits == 1 and rev_hits == 1:
                # 5. Apply Tailing Logic (if specific)
                fwd_rc = str(Seq(fwd_seq).reverse_complement())
                rev_rc = str(Seq(rev_seq).reverse_complement())
                pair_rank_str = f"{i} (Strategy {strategy_index})"

                csv_row = {
                    'target_id': target_id, 'pair_rank': pair_rank_str, 
                    'fwd_primer_seq': fwd_seq, 'rev_primer_seq': rev_seq,
                    'fwd_primer_tailed': fwd_rc + FWD_TAIL,
                    'rev_primer_tailed': rev_rc + REV_TAIL,
                    'fwd_primer_tm': f"{primer_results[f'PRIMER_LEFT_{i}_TM']:.2f}",
                    'rev_primer_tm': f"{primer_results[f'PRIMER_RIGHT_{i}_TM']:.2f}",
                    'amplicon_size': primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                    'specificity_hits': f"F:{fwd_hits}, R:{rev_hits}"
                }
                
                product_size = primer_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                fwd_start_in_gene = primer_results[f'PRIMER_LEFT_{i}'][0]
                amp_start = target_info['start'] - 1 + fwd_start_in_gene
                amp_end = amp_start + product_size
                bed_row = f"{target_info['contig']}\t{amp_start}\t{amp_end}\t{target_id}_amplicon_{pair_rank_str}\t0\t+\n"
                
                all_specific_pairs_for_target.append({'csv_row': csv_row, 'bed_row': bed_row})

    if not all_specific_pairs_for_target:
        return (target_id, [], f"Could not find a specific primer pair for {target_id} after all strategies.")
    
    return (target_id, all_specific_pairs_for_target, None) # Return target_id, list of results, and no error

def _check_pair_compatibility(p1_name, p1_seq, p2_name, p2_seq):
    """Helper function to check a single pair for compatibility."""
    # Define a single, conservative Delta G threshold (in kcal/mol).
    DG_THRESHOLD = -6.0
    
    try:
        # Check Heterodimer (cross-dimer)
        result = primer3.calc_heterodimer(p1_seq, p2_seq)
        # Fix: divide by 1000.0 to get kcal/mol
        dg_val = result.dg / 1000.0
        if dg_val < DG_THRESHOLD:
            return (
                f"Potential cross-dimer between {p1_name} and {p2_name} "
                f"(dG: {dg_val:.2f})"
            )
    except Exception as e:
        print(f"Warning: Could not calculate heterodimer for {p1_name}/{p2_name}. Error: {e}")
    
    return None

def find_compatible_set(target_to_primers_map, max_iterations=50):
    """
    Attempts to find a compatible set of primers by iterating and
    swapping out clashing pairs.
    
    NEW (v3.2): If a perfect set isn't found, it returns the "best available"
    imperfect set (the one with the fewest clashes).
    """
    print("\n--- Starting Multiplex Compatibility Check & Auto-Healing ---")
    
    if len(target_to_primers_map) < 2:
        print("Not enough targets for a multiplex check. Using best primers.")
        # Just return the best primer for each target
        final_set = [primers[0] for primers in target_to_primers_map.values() if primers]
        return final_set, []

    # 1. Initialize: Start with the "best" primer (index 0) for every target
    current_indices = {target: 0 for target in target_to_primers_map}
    
    # --- NEW v3.2 Logic ---
    # Keep track of the best set found so far
    best_set_indices = current_indices.copy()
    min_clashes_found = float('inf')
    best_clash_list = []
    # --- End NEW Logic ---
    
    for iteration in range(max_iterations):
        current_primer_set = [] # List of (primer_name, seq, target_id)
        
        # 2. Build the current set based on current_indices
        valid_set = True
        for target_id, index in current_indices.items():
            if index >= len(target_to_primers_map[target_id]):
                print(f"Error: Ran out of primer alternatives for {target_id} during iteration {iteration+1}.")
                valid_set = False
                break
            
            primer_data = target_to_primers_map[target_id][index]
            csv_row = primer_data['csv_row']
            fwd_name = f"{csv_row['target_id']}_F_{csv_row['pair_rank']}"
            rev_name = f"{csv_row['target_id']}_R_{csv_row['pair_rank']}"
            
            current_primer_set.append((fwd_name, csv_row['fwd_primer_seq'], target_id))
            current_primer_set.append((rev_name, csv_row['rev_primer_seq'], target_id))
        
        if not valid_set:
            # We ran out of options for a target, so we can't continue this path.
            # We must return the best set we've seen so far.
            print(f"Stopping search. Returning best set found (with {min_clashes_found} clashes).")
            final_set = [target_to_primers_map[t][i] for t, i in best_set_indices.items() if i < len(target_to_primers_map[t])]
            return final_set, best_clash_list

        # 3. Check for clashes in the current set
        clashes = {} 
        all_incompatible_pairs = []

        # Check for cross-dimers
        for (p1_name, p1_seq, p1_target), (p2_name, p2_seq, p2_target) in itertools.combinations(current_primer_set, 2):
            if p1_target == p2_target:
                continue
            error_msg = _check_pair_compatibility(p1_name, p1_seq, p2_name, p2_seq)
            if error_msg:
                all_incompatible_pairs.append(error_msg)
                clashes[p1_target] = clashes.get(p1_target, 0) + 1
                clashes[p2_target] = clashes.get(p2_target, 0) + 1

        # 4. Check for Self-Dimers
        for name, seq, target_id in current_primer_set:
             try:
                result = primer3.calc_homodimer(seq)
                dg_val = result.dg / 1000.0
                if dg_val < -6.0:
                    error_msg = f"Potential self-dimer in {name} (dG: {dg_val:.2f})"
                    all_incompatible_pairs.append(error_msg)
                    clashes[target_id] = clashes.get(target_id, 0) + 1
             except Exception as e:
                print(f"Warning: Could not calculate homodimer for {name}. Error: {e}")

        # --- NEW v3.2 Logic ---
        # 5. Evaluate and 'Heal' or Exit
        num_clashes = len(all_incompatible_pairs)
        
        if num_clashes < min_clashes_found:
            min_clashes_found = num_clashes
            best_set_indices = current_indices.copy()
            best_clash_list = all_incompatible_pairs
        # --- End NEW Logic ---

        if not clashes:
            print(f"SUCCESS: Found a compatible set in {iteration + 1} iterations.")
            final_set = [target_to_primers_map[t][i] for t, i in current_indices.items()]
            return final_set, []
        
        worst_target = max(clashes, key=clashes.get)
        print(f"  -> Iteration {iteration+1}: Found {num_clashes} clashes. Worst offender: {worst_target} ({clashes[worst_target]} clashes).")
        
        current_indices[worst_target] += 1
        
    # 6. If loop finishes, we failed. Return the BEST set we found.
    print(f"Failed to find a perfect set after {max_iterations} iterations.")
    print(f"Returning best available set with {min_clashes_found} potential clashes.")
    final_set = [target_to_primers_map[t][i] for t, i in best_set_indices.items() if i < len(target_to_primers_map[t])]
    return final_set, best_clash_list

def run_design_mode(args):
    """Runs the script in 'Full Design' mode."""
    print("Running in 'Full Design' mode...")
    failed_targets_initial = []

    try:
        # 1. Load shared data ONCE
        create_blast_db_if_needed(args.genome, args.blast_db)
        cleaned_genome_path = f"{os.path.splitext(args.genome)[0]}.cleaned.fna"
        print(f"Loading CLEANED genome from '{cleaned_genome_path}'...")
        genome_records = SeqIO.to_dict(SeqIO.parse(cleaned_genome_path, "fasta"))
        print(f"Parsing GFF file '{args.gff}'...")
        gene_coords = parse_gff(args.gff)
        print(f"Reading target IDs from '{args.target_file}'...")
        target_ids = read_lines_from_file(args.target_file)

        # 2. Create a "partial" function for parallel processing
        process_func = partial(process_single_target,
                               genome_records=genome_records,
                               gene_coords=gene_coords,
                               blast_db=args.blast_db)

        # 3. Create a process pool and run the tasks in parallel
        num_workers = os.cpu_count()
        print(f"Starting parallel processing with {num_workers} workers for {len(target_ids)} targets...")
        
        results = []
        with multiprocessing.Pool(processes=num_workers) as pool:
            results = list(tqdm(pool.imap_unordered(process_func, target_ids), total=len(target_ids), desc="Designing Primers"))

        # 4. Collect results into a map: {target_id -> [list_of_primer_options]}
        target_to_primers_map = {}
        for target_id, primer_list, error in results:
            if primer_list:
                target_to_primers_map[target_id] = primer_list
            else:
                if error: 
                    failed_targets_initial.append(error)
        
        if not target_to_primers_map:
             print("\nNo specific primers were found for any target. Exiting.")
             return

        # 5. Run the Auto-Healing Resolver
        final_compatible_set, final_warnings = find_compatible_set(target_to_primers_map)

        if final_warnings:
            print("\n--- Final Compatibility Warnings ---")
            for warning in sorted(list(set(final_warnings))): # Use set to remove duplicates
                print(f"  - {warning}")
        
        if failed_targets_initial:
            print("\n--- Targets That Failed Initial Design ---")
            for reason in sorted(list(set(failed_targets_initial))): # Use set to remove duplicates
                print(f"  - {reason}")

        # 6. Write final files
        if final_compatible_set:
            final_csv_rows = [row['csv_row'] for row in final_compatible_set]
            final_bed_rows = [row['bed_row'] for row in final_compatible_set]
            write_design_output_files(final_csv_rows, final_bed_rows, args.output_prefix)
        else:
            print("\nCould not determine a final compatible set. No files will be written.")

    except (FileNotFoundError, subprocess.CalledProcessError, RuntimeError, ValueError) as e:
        print(f"\nAn error occurred: {e}")
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")

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
        
        # Logic fix: Primer first, then tail
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

# --- Main Execution Block (Router) ---

def main():
    parser = argparse.ArgumentParser(description="Design specific, adapter-tailed primers for multiple targets.")
    
    # Mode 1: Full Design
    design_group = parser.add_argument_group('Mode 1: Full Design Pipeline')
    design_group.add_argument('--genome', help="Path to the reference genome in FASTA format.")
    design_group.add_argument('--gff', help="Path to a GFF file for gene coordinate lookups.")
    design_group.add_argument('--target-file', help="Path to a text file with one target gene ID per line.")
    design_group.add_argument('--blast-db', help="Prefix for the local BLAST database.")
    
    # Mode 2: Tail-Only
    tail_group = parser.add_argument_group('Mode 2: Tail-Only Utility')
    tail_group.add_argument('--tail-fwd-file', help="Path to a text file with one forward primer per line.")
    tail_group.add_argument('--tail-rev-file', help="Path to a text file with one reverse primer per line.")
    
    # Shared arguments
    parser.add_argument('--output-prefix', default='final_primers', help="Prefix for output files.")
    
    args = parser.parse_args()

    # --- Route to the correct mode ---
    if args.genome and args.gff and args.target_file and args.blast_db:
        run_design_mode(args)
    elif args.tail_fwd_file and args.tail_rev_file:
        run_tail_only_mode(args)
    else:
        print("Error: You must provide the correct arguments for a mode.")
        print("\nFor 'Full Design' mode, you MUST provide:")
        print("  --genome, --gff, --target-file, and --blast-db")
        print("\nFor 'Tail-Only' mode, you MUST provide:")
        print("  --tail-fwd-file and --tail-rev-file")
        parser.print_help()

if __name__ == "__main__":
    main()

