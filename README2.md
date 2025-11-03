# Multiplex PCR Primer Designer (v4.0)

This is a powerful command-line tool for designing robust, multiplex-ready PCR primer panels for amplicon-based next-generation sequencing.

It is an end-to-end solution that takes a genome, a list of gene targets, and produces a final, validated, and sequencing-ready primer list. It is built to handle the major challenges of multiplex PCR, including primer specificity, primer-dimer interactions, hairpin formation, and on-target cross-priming (tiling).

### Core Features
* **Two-Mode Operation:** Functions as a full "Design" pipeline or a simple "Tail-Only" utility.
* **Parallel Processing:** Uses all available CPU cores to design large panels quickly.
* **Multi-Check Validation:** All primers are validated for:
    1.  **Genome Specificity** (via local BLAST)
    2.  **Hairpin Structures** (via `primer3.calc_hairpin`)
    3.  **Dimer-Compatibility** (via `primer3.calc_homodimer/calc_heterodimer`)
* **Smart Retry Logic:** Automatically retries failed targets with different design parameters.
* **"Auto-Healing" Resolver:** Intelligently finds the "best-available" compatible primer set for a single pool, even if a perfect set is impossible.
* **Tiled Design Mode:** Automatically solves on-target cross-priming by sorting tiled amplicons into separate, interleaved pools.
* **Actionable Output:** Generates a `.csv` for your order, a `.bed` for visualization, and a `.log.txt` for full traceability of design warnings.

---

## Files in This Repository

* `design.py`: The main, executable Python script.
* `requirements.txt`: A list of all required Python libraries.
* `.gitignore`: A list of files (like results and BLAST databases) to ignore in version control.
* `README.md`: This file.

---

## Installation

**1. Clone the Repository**
```bash
git clone [https://github.com/YOUR_USERNAME/primer_designer.git](https://github.com/YOUR_USERNAME/primer_designer.git)
cd primer_designer
```

**2. Install Python Dependencies**
This tool requires `biopython`, `primer3-py`, and `tqdm`.
```bash
pip install -r requirements.txt
```

**3. Install NCBI BLAST+**
This tool requires the BLAST+ command-line tools to be installed on your system.
* **Download:** [NCBI BLAST+ Download Page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (Download the `...-win64.exe` installer for Windows).
* **Install:** Run the installer. **Crucially, ensure you check the box "Add BLAST+ to the system PATH."**
* **Verify:** Open a *new* terminal/PowerShell and type `blastn -version`. You should see a version number.

---

## Usage

The script operates in two main modes, which it selects based on the arguments you provide.

### Mode 1: Full Design Pipeline

This is the main mode. It designs, validates, and tails primers from scratch. It requires you to specify a `--design-type`.

**Required Arguments:**
* `--genome`: Path to your reference genome (FASTA format).
* `--gff`: Path to your gene annotation file (GFF format).
* `--target-file`: Path to a simple `.txt` file with one gene name per line.
* `--blast-db`: A prefix for your BLAST database (e.g., "ecoli_db"). The script will automatically build this database from your genome if it doesn't exist.
* `--design-type`: `sparse` or `tiled`. This is the most important setting.

#### Design Type: `sparse` (Single-Pool)
This mode attempts to create the best possible primer panel for a **single PCR tube**. It uses the "Auto-Healing" logic to find the set with the *fewest possible clashes*.

* **Best For:** Convenience, non-tiled targets, and lower-plex panels (< 20 targets).
* **Output:** A single set of files (`.csv`, `.bed`, `.log.txt`). It will warn you if the best set is still high-risk.

**Example Command:**
```powershell
python .\design.py --genome "ecoli_genome.fna" --gff "genomic.gff" --target-file "target_genes.txt" --blast-db "ecoli_db" --design-type sparse --output-prefix "ecoli_10plex_sparse_run"
```

#### Design Type: `tiled` (Multi-Pool)
This mode is designed for **high-plex panels** or **tiled amplicons** where on-target cross-priming is a risk. It automatically solves this by splitting targets into two (or more) interleaved pools.

* **Best For:** High-plexity panels (> 20 targets), tiled amplicons, and guaranteeing a 0-clash, high-performance result.
* **Output:** Two sets of files: `_pool_A.csv`, `_pool_A.bed`, etc., and `_pool_B.csv`, `_pool_B.bed`, etc.

**Example Command:**
```powershell
python .\design.py --genome "ecoli_genome.fna" --gff "genomic.gff" --target-file "target_genes.txt" --blast-db "ecoli_db" --design-type tiled --output-prefix "ecoli_10plex_tiled_run"
```

---

### Mode 2: Tail-Only Utility

This mode bypasses all design and validation. It simply takes existing primer lists and adds your hardcoded sequencing tails (using the reverse-complement logic).

**Example Command:**
```powershell
python .\design.py --tail-fwd-file "my_fwd_primers.txt" --tail-rev-file "my_rev_primers.txt" --output-prefix "my_tailed_primers"

