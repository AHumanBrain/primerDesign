# Primer Design and Tailing Tool (v6.0)

This Python tool provides a robust pipeline for designing specific, adapter-tailed PCR primers for multiple target genes against a reference genome. It also includes a utility mode to simply add adapter tails to existing primer sequences. The core design logic leverages `primer3-py` for primer design and NCBI BLAST+ for specificity checks, with advanced features for multiplex compatibility assessment and automated pool generation.

## Table of Contents
- [Core Features](#core-features)
- [Included Files](#included-files)
- [Installation](#installation)
- [Workflow Logic (Design Mode)](#workflow-logic-design-mode)
- [Usage](#usage)
  - [Mode 1: Full Design Pipeline](#mode-1-full-design-pipeline)
  - [Mode 2: Tail-Only Utility](#mode-2-tail-only-utility)
- [Understanding the Output](#understanding-the-output)
- [Design Constants & Customization](#design-constants--customization)

## Core Features

* **Targeted Primer Design:** Designs primers for specified gene IDs from a GFF file, using sequences extracted from a reference FASTA.
* **BLAST Specificity Check:** Ensures primers are specific to their intended target in the genome by requiring single, 100% identity hits (can be relaxed with `--force-multiprime`).
* **3' End Stability Filtering:** Filters out primers prone to self-dimerization, hairpin formation, or cross-dimerization at the 3' end, a critical step for multiplex PCR.
* **Multi-Strategy Design:** Employs several `primer3` parameter sets (e.g., different amplicon sizes, relaxed Tm ranges) to find suitable primers for challenging targets, prioritizing "ideal" conditions first.
* **Adapter Tailing:** Automatically adds predefined adapter sequences to the 5' end of the reverse-complement of designed primers.
* **Multiplex Compatibility (`find_compatible_set`):**
    * **Auto-Healing Algorithm:** Iteratively selects primer pairs to minimize inter-primer 3' end cross-dimerization, aiming for zero clashes.
    * **"Best Available" Set:** If a perfect set isn't found, it returns the set with the fewest and weakest clashes within a configurable number of iterations.
* **Smart Pooling Strategy:**
    * **Overlap Detection:** Automatically detects if amplicons for a set of targets overlap (indicating a "tiled" design).
    * **Dynamic Pooling:** For tiled designs or by default, it finds the *minimum number of compatible pools* needed for all targets.
    * **Forced Single Pool:** Option to force a single multiplex pool, even if clashes are present, returning the best possible compromise.
* **Parallel Processing:** Utilizes `multiprocessing` to speed up primer design for many targets.
* **Tail-Only Mode:** A convenient utility to add standard adapter tails to a list of pre-existing forward and reverse primers.
* **Comprehensive Output:** Generates CSV files with primer details (sequences, Tms, flags, specificity), BED files for amplicon visualization, and a detailed log file with warnings and failed targets.

## Included Files

* `design_v6.py`: The main Python script containing all functionalities.

## Installation

This tool requires Python 3.7+ and several bioinformatics libraries and tools.

1.  **Clone the Repository (or download `design_v6.py`):**
    ```bash
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    cd your-repo-name
    ```

2.  **Install Python Dependencies:**
    It's highly recommended to use a virtual environment.
    ```bash
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    pip install biopython primer3-py tqdm
    ```

3.  **Install NCBI BLAST+:**
    Download and install the appropriate NCBI BLAST+ executables for your operating system from the [NCBI FTP site](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
    Ensure that `makeblastdb` and `blastn` are in your system's PATH. You can test this by opening a new terminal and typing `makeblastdb -help`.

    * **Linux/macOS:** You might install via a package manager (`sudo apt-get install ncbi-blast+` or `brew install blast`).
    * **Windows:** Add the installation directory (e.g., `C:\ncbi-blast-2.13.0+\bin`) to your system's PATH environment variable.

## Workflow Logic (Design Mode)

The design pipeline follows a structured approach to ensure high-quality, multiplex-compatible primers:

```mermaid
graph TD
    A[Start Design Mode] --> B{Load Genome, GFF, Targets};
    B --> C{Create/Check BLAST DB};
    C --> D{Parallel Primer Design for Each Target};
    D --> D0[Strategy 0: Ideal Tm (59-61C), Prod Size (150-250bp)];
    D --> D1[Strategy 1: Longer Prod Size (250-350bp)];
    D --> D2[Strategy 2: Shorter Prod Size (100-150bp)];
    D --> D3[Strategy 3: Relaxed Tm (55-65C)];
    D0 & D1 & D2 & D3 --> E{Specificity & 3' End Stability Check};
    E --> F{Collect All Valid Primer Options per Target};
    F --> G{Check for Amplicon Overlaps (Best Primary Set)};
    G -- Overlaps Detected (Tiled Design) --> H1[Force Multi-Pool Design];
    G -- No Overlaps (Sparse Design) --> H2{User `--force-single-pool`?};
    H2 -- No (Default) --> H1;
    H2 -- Yes --> I[Attempt Single Pool Compatibility];
    H1 --> J{Iteratively Find Minimum N Compatible Pools};
    J --> K[Output N Pool Files];
    I --> L[Output Single Pool File];
    K & L --> M[Generate Log & Warnings];
    M --> N[End];
    E -- No Valid Pair --> O[Log Failed Target (Initial Design)];
    O --> F;
```

## Usage

Run the script from your terminal. Use `python design_v6.py --help` for a full list of arguments.

### Mode 1: Full Design Pipeline

This mode requires a genome FASTA, GFF annotation, a list of target gene IDs, and a BLAST database prefix.

**Required Arguments:**

* `--genome <path/to/genome.fasta>`: Path to the reference genome FASTA file.
* `--gff <path/to/annotation.gff>`: Path to the GFF file for gene coordinates.
* `--target-file <path/to/targets.txt>`: Path to a text file, one target gene ID per line.
* `--blast-db <prefix>`: Prefix for the local BLAST database (e.g., `ecoli_db`). If the database doesn't exist, it will be created.

**Optional Arguments:**

* `--output-prefix <prefix>`: Prefix for output files (default: `final_primers`).
* `--force-single-pool`: If specified, the script will attempt to design for a single multiplex pool, even if compatibility issues exist. By default, it finds the minimum number of compatible pools.
* `--force-multiprime`: Allows primers that hit multiple identical locations in the genome to be considered valid. Useful for multi-copy genes (e.g., 16S rRNA) where a single specific hit isn't expected. Default behavior is strict single-hit specificity.
* `--max-compatibility-iterations <int>`: Maximum attempts for the auto-healing algorithm to find a compatible set (default: 100). Increase for very complex panels if clashes persist.

**Example Command (Basic):**

```bash
python design_v6.py \
    --genome ./data/ecoli.fasta \
    --gff ./data/ecoli.gff \
    --target-file ./data/target_genes.txt \
    --blast-db ecoli_blast_db \
    --output-prefix my_panel_design
```

**Example Command (Forced Single Pool, Multi-Copy Targets):**

```bash
python design_v6.py \
    --genome ./data/ecoli.fasta \
    --gff ./data/ecoli.gff \
    --target-file ./data/16s_rRNA_targets.txt \
    --blast-db ecoli_blast_db \
    --force-single-pool \
    --force-multiprime \
    --max-compatibility-iterations 200 \
    --output-prefix 16s_multiplex
```

### Mode 2: Tail-Only Utility

This mode allows you to add adapter tails to existing primer sequences.

**Required Arguments:**

* `--tail-fwd-file <path/to/fwd_primers.txt>`: Path to a text file with one forward primer sequence per line.
* `--tail-rev-file <path/to/rev_primers.txt>`: Path to a text file with one reverse primer sequence per line.

**Optional Arguments:**

* `--output-prefix <prefix>`: Prefix for output files (default: `final_primers`).

**Example Command:**

```bash
python design_v6.py \
    --tail-fwd-file ./my_fwd_primers.txt \
    --tail-rev-file ./my_rev_primers.txt \
    --output-prefix tailed_my_primers
```

## Understanding the Output

For **Design Mode**, the script will generate the following files (using `my_panel_design` as `output-prefix`):

* **`my_panel_design.csv` (or `my_panel_design_pool_N.csv`):**
    A CSV file containing detailed information for each designed primer pair:

    * `target_id`: The ID of the gene for which the primers were designed.
    * `pair_rank`: The rank of the primer pair from Primer3, plus the strategy used (e.g., `0 (Strategy 0)`).
    * `flags`: Any warnings or issues (e.g., `Low_Tm`, `High_Tm` if outside the `IDEAL_TM_MIN/MAX` range, `Multi_Hit_FWD:X_REV:Y` if `--force-multiprime` was used). `OK` indicates no flags.
    * `fwd_primer_tailed`: The 5' reverse-complement of the forward primer + FWD adapter tail.
    * `rev_primer_tailed`: The 5' reverse-complement of the reverse primer + REV adapter tail.
    * `fwd_primer_seq`: The original forward primer sequence as designed by Primer3.
    * `rev_primer_seq`: The original reverse primer sequence as designed by Primer3.
    * `fwd_primer_tm`: Melting temperature of the forward primer.
    * `rev_primer_tm`: Melting temperature of the reverse primer.
    * `amplicon_size`: The expected size of the PCR product.
    * `specificity_hits`: Number of 100% identical hits for FWD and REV primers in the genome (e.g., `F:1, R:1`).

* **`my_panel_design.bed` (or `my_panel_design_pool_N.bed`):**
    A BED file listing the genomic coordinates of each designed amplicon. This can be loaded into genome browsers (e.g., IGV, UCSC Genome Browser) to visualize amplicon locations.

* **`my_panel_design.log.txt`:**
    A comprehensive log file containing:

    * Details of the design process.
    * Any final compatibility warnings for the chosen primer set (from `find_compatible_set`).
    * A list of targets for which no suitable primer pairs could be designed, along with the reasons.

For **Tail-Only Mode**, the script will generate:

* **`tailed_my_primers.csv`:**
    A CSV file with the original and tailed primer sequences:
    * `pair_id`: An incremental ID for each primer pair.
    * `fwd_primer_tailed`: The 5' reverse-complement of the forward primer + FWD adapter tail.
    * `rev_primer_tailed`: The 5' reverse-complement of the reverse primer + REV adapter tail.
    * `fwd_primer_seq`: The original forward primer sequence.
    * `rev_primer_seq`: The original reverse primer sequence.

## Design Constants & Customization

The following constants are defined at the top of `design_v6.py` and can be adjusted directly in the script to fine-tune primer design behavior and criteria:

```python
# --- Hardcoded Adapter Tails ---
FWD_TAIL = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
REV_TAIL = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

# --- Design Constants ---
END_STABILITY_DG_THRESHOLD = -9.0  # (kcal/mol) Max allowed 3' end interaction strength. More negative is worse.
IDEAL_TM_MIN = 59.0                # Minimum Tm for 'OK' flag; ALSO used as PRIMER_MIN_TM for Strategy 0
IDEAL_TM_MAX = 61.0                # Maximum Tm for 'OK' flag; ALSO used as PRIMER_MAX_TM for Strategy 0
MAX_CLASH_RECOMMENDATION = 5       # Number of allowed clashes before a warning is issued for single pools.

# Primer3 global arguments are set within design_primers_for_sequence.
# These define the initial target parameters for Primer3.
# Important parameters here (refer to Primer3 documentation for full list):
# 'PRIMER_OPT_SIZE': 20,         # Target primer length
# 'PRIMER_MIN_SIZE': 19,         # Minimum primer length
# 'PRIMER_MAX_SIZE': 22,         # Maximum primer length
# 'PRIMER_OPT_TM': 62.5,         # Primer3 actively tries to design primers with this Tm
# 'PRIMER_MIN_GC': 40.0,         # Minimum GC content percentage
# 'PRIMER_MAX_GC': 60.0,         # Maximum GC content percentage
# 'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], # Default amplicon size range
# 'PRIMER_NUM_RETURN': 20        # Number of primer pairs Primer3 attempts to return per sequence
```

The `strategies` list, defined within the `process_single_target` function, provides sequential alternative parameter sets for Primer3 to try if the initial attempt fails. Modifying these allows you to define a precise retry logic:

```python
# --- Inside process_single_target function ---
    # 1. Define Retry Strategies
    strategies = [
        # Strategy 0: Your primary, tightest attempt, using the IDEAL Tms
        {'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_MIN_TM': IDEAL_TM_MIN, 'PRIMER_MAX_TM': IDEAL_TM_MAX}, 
        # Strategy 1: Longer product size range, with slightly relaxed Tms from the ideal
        {'PRIMER_PRODUCT_SIZE_RANGE': [[250, 350]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}, 
        # Strategy 2: Shorter product size range, with slightly relaxed Tms from the ideal
        {'PRIMER_PRODUCT_SIZE_RANGE': [[100, 150]], 'PRIMER_MIN_TM': 57.0, 'PRIMER_MAX_TM': 63.0}, 
        # Strategy 3: Default product size, but with the most relaxed Tm range
        {'PRIMER_PRODUCT_SIZE_RANGE': [[150, 250]], 'PRIMER_MIN_TM': 55.0, 'PRIMER_MAX_TM': 65.0}, 
    ]
```

**Understanding TM Settings:**
* **`PRIMER_OPT_TM` (in `global_args`):** This is the single temperature Primer3 *tries hardest* to achieve. For robust multiplexing, aiming for a slightly higher optimum (e.g., 62.5°C) within an acceptable range (e.g., 59-65°C) can be beneficial.
* **`PRIMER_MIN_TM` / `PRIMER_MAX_TM` (in `global_args` and `strategies`):** These define the *hard boundaries* for Primer3's search for Tm. A primer designed outside these will be discarded by Primer3.
* **`IDEAL_TM_MIN` / `IDEAL_TM_MAX` (top-level constants):** These serve a dual purpose. They are used as the `PRIMER_MIN_TM` / `PRIMER_MAX_TM` for the initial, most stringent `Strategy 0`. Additionally, they define the bounds for post-design flagging (`Low_Tm` / `High_Tm`). This means a primer designed by a more relaxed strategy (e.g., `57.0-63.0`) might still receive a `High_Tm` flag if its Tm is above `IDEAL_TM_MAX` (e.g., 61.0°C), indicating it's acceptable by relaxed criteria but not "ideal."
````http://googleusercontent.com/image_generation_content/9
