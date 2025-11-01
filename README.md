# Multiplex PCR Primer Designer

This is a command-line tool for designing multiplex-ready PCR primers. It provides a complete end-to-end workflow, from target selection to a final, validated set of primers with sequencing adapter tails.

The script is built to be robust, fast, and specific, using a combination of Primer3 for design, NCBI BLAST+ for specificity checking, and an "auto-healing" algorithm to resolve primer-dimer incompatibilities.

## Features

* **Two Modes of Operation:**
    1.  **Full Design Mode:** Designs, validates, and tails new primers from scratch for a list of gene targets.
    2.  **Tail-Only Mode:** Bypasses design/validation and simply adds sequencing tails to a pre-designed list of primers.
* **Intelligent Design Pipeline:**
    * **Parallel Processing:** Uses all available CPU cores to process large primer panels in parallel.
    * **Specificity Guarantee:** Uses a local NCBI BLAST database to ensure every primer pair binds *only* to its intended target.
    * **Smart Retry:** Automatically tries multiple design strategies (e.g., different amplicon sizes, Tms) for challenging targets that fail on the first pass.
* **Advanced Multiplex Validation:**
    * **Cross-Compatibility Check:** Checks all primers in the final set against each other for potential primer-dimers ($\Delta G$).
    * **Auto-Healing Algorithm:** If clashes are found, the script intelligently swaps out "problem" primers with the next-best specific alternatives and re-checks, attempting to find a fully compatible set.
    * **Best-Available Set:** If a 100% perfect set can't be found, the script outputs the "best available" set (the one with the fewest clashes) along with a clear warning list.
* **Automated Tailing:** Applies custom sequencing tails using reverse-complement logic, ready for ordering.
* **Clear Outputs:** Generates a `.csv` file for analysis and a `.bed` file for visualizing amplicons in a genome browser like IGV.

## Requirements

1.  **Python 3.7+**
2.  **NCBI BLAST+** (must be installed and in your system's `PATH`)
3.  **Python Libraries** (install with `pip install -r requirements.txt`):
    * `biopython`
    * `primer3-py`
    * `tqdm`

## `requirements.txt`

```
biopython
primer3-py
tqdm
```

## Initial Setup: Building the BLAST Database

Before you can run the design pipeline, you must create a local BLAST database from your genome. The script can do this for you, but it's best to run it once to be sure.

1.  Place your genome `your_genome.fna` in your project folder.
2.  Run the `makeblastdb` command (part of BLAST+):

```powershell
makeblastdb -in "your_genome.fna" -dbtype nucl -out "your_blast_db_name" -parse_seqids
```

*Note: The script also has a built-in function (`create_blast_db_if_needed`) that will try to do this automatically if it can't find the database, including cleaning FASTA headers that would otherwise cause `makeblastdb` to fail.*

## Usage

The script operates in two distinct modes.

---

### Mode 1: Full Design Pipeline

This is the main workflow. You provide a genome, gene annotations, and a list of target genes, and the script generates a final, multiplex-compatible primer set.

**Example Command:**

```powershell
python .\design.py --genome "ecoli_genome.fna" `
                  --gff "genomic.gff" `
                  --target-file "target_genes.txt" `
                  --blast-db "ecoli_db" `
                  --output-prefix "ecoli_v3_final_run"
```

**Example `target_genes.txt`:**

```
gyrA
recA
dnaA
rpsL
trpA
ampC
lacZ
rpoB
fliC
dnaB
```

**Example Output Log:**

```
Running in 'Full Design' mode...
BLAST database 'ecoli_db' already exists. Skipping creation.
Loading CLEANED genome from 'ecoli_genome.cleaned.fna'...
Parsing GFF file 'genomic.gff'...
Reading target IDs from 'target_genes.txt'...
Starting parallel processing with 8 workers for 10 targets...
Designing Primers: 100%|██████████| 10/10 [00:01<00:00, 6.41it/s]

--- Starting Multiplex Compatibility Check & Auto-Healing ---
  -> Iteration 1: Found 22 clashes. Worst offender: rpoB (7 clashes).
  -> Iteration 2: Found 19 clashes. Worst offender: recA (6 clashes).
  ...
  -> Iteration 49: Found 8 clashes. Worst offender: dnaA (3 clashes).
  -> Iteration 50: Found 11 clashes. Worst offender: dnaA (6 clashes).
Failed to find a perfect set after 50 iterations.
Returning best available set with 8 potential clashes.

--- Final Compatibility Warnings ---
  - Potential cross-dimer between ampC_R_1 (Strategy 0) and recA_R_2 (Strategy 0) (dG: -6.21)
  - Potential cross-dimer between dnaA_F_0 (Strategy 0) and trpA_R_0 (Strategy 0) (dG: -6.21)
  ... (6 more warnings) ...

Results for 10 specific targets saved to 'ecoli_v3.1_autoheal_run.csv' and 'ecoli_v3.1_autoheal_run.bed'
```

---

### Mode 2: Tail-Only Utility

Use this mode if you *already* have primer sequences and simply want to apply the hardcoded adapter tails. It bypasses all design, specificity, and compatibility checks.

**Example Command:**

```powershell
python .\design.py --tail-fwd-file "my_fwd_primers.txt" `
                  --tail-rev-file "my_rev_primers.txt" `
                  --output-prefix "tailed_primer_set"
```

**Example `my_fwd_primers.txt`:**

```
GCTATCACCCAGTTTGATCG
AACTCACTTCGGTCAGGTCG
...
```

---

### Understanding the Output CSV

The main output file (`.csv`) contains the final primer set. Here are the key columns:

* **`target_id`**: The gene name you provided.
* **`pair_rank`**: The rank from Primer3. `0 (Strategy 0)` is the "best" primer pair from the default strategy. `1 (Strategy 2)` would be the *second-best* primer pair from the *third* (index 2) retry strategy.
* **`fwd_primer_seq`**: The specific forward primer sequence.
* **`rev_primer_seq`**: The specific reverse primer sequence.
* **`fwd_primer_tailed`**: The final, ready-to-order sequence for the forward primer.
* **`rev_primer_tailed`**: The final, ready-to-order sequence for the reverse primer.
* **`fwd_primer_tm`**: Melting temperature of the specific forward primer.
* **`rev_primer_tm`**: Melting temperature of the specific reverse primer.
* **`amplicon_size`**: The length of the amplicon *before* tails are added.
* **`specificity_hits`**: The result of the BLAST check. **`F:1, R:1`** is the ideal, meaning the forward primer and reverse primer each bind exactly once in the genome.
