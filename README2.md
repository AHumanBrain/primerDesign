<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>README Copier (v4.3 - FINAL FIX 2)</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <style>
        body { font-family: 'Inter', sans-serif; }
        /* Ensure textarea respects markdown formatting for copying */
        textarea { 
            white-space: pre; 
            overflow-wrap: normal; 
            overflow-x: scroll; 
            font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
            font-size: 0.875rem;
            line-height: 1.25rem;
        }
    </style>
</head>
<body class="bg-slate-100 p-4 md:p-8">
    <div class="max-w-4xl mx-auto bg-white p-6 md:p-8 rounded-lg shadow-lg">
        <h1 class="text-2xl font-bold text-slate-800">Updated README.md (v4.3 - FINAL FIX 2)</h1>
        <p class="text-slate-600 mt-2">The text box below contains the **complete and corrected** raw Markdown for your project's README. Use the button to copy the *entire* contents to your clipboard, then paste it directly into your `README.md` file on GitHub.</p>
        
        <textarea id="readmeText" class="w-full h-96 mt-4 p-3 border border-slate-300 rounded-md bg-slate-50" readonly>
# Multiplex PCR Primer Designer (v4.3)

This is a powerful, fully-automated command-line tool for designing robust, multiplex-ready PCR primer panels for amplicon-based next-generation sequencing.

It is an end-to-end solution that takes a genome, a list of gene targets, and produces a final, validated, and sequencing-ready primer list. The script automatically detects whether your design is "sparse" or "tiled" (overlapping) and applies the optimal design strategy to ensure the highest probability of success at the lab bench.

### Core Features
* **Fully-Automated Design Logic:** Automatically detects tiled/overlapping amplicons and forces a multi-pool design to prevent on-target cross-priming.
* **Two-Mode Operation:** Functions as a full "Design" pipeline or a simple "Tail-Only" utility.
* **Parallel Processing:** Uses all available CPU cores to design large panels quickly.
* **Multi-Check Validation:** All primers are validated for:
    1.  **Genome Specificity** (via local BLAST)
    2.  **Hairpin Structures** (via `primer3.calc_hairpin`)
    3.  **Dimer-Compatibility** (via `primer3.calc_homodimer/calc_heterodimer`)
* **Smart Retry Logic:** Automatically retries failed targets with different design parameters (size, $T_m$, etc.).
* **Intelligent Auto-Pooling (v4.2):** For sparse or tiled designs, this is the default mode. It finds the *minimum number of compatible pools* (N-pools) required to achieve a **perfect, 0-clash** design, minimizing your lab burden.
* **Forced Single-Pool Mode (v3.2):** For convenience, you can force the script to find the "best-available" (but possibly imperfect) set for a single tube.
* **Actionable Output:** Generates a `.csv` for your order (with optimization flags), a `.bed` for visualization, and a `.log.txt` for full traceability.

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

## Workflow Logic (v4.3)

The script's automated decision-making follows this path:

```mermaid
graph TD
    A[Start: Run design.py] --> B[Run Parallel Design v2.1]
    B --> C["Find all specific, hairpin-free<br>primer options for all targets"]
    C --> D{Check if 'best' amplicons overlap}
    D -- Yes, Tiled --> E[Force Auto-N-Pool Logic]
    D -- No, Sparse --> F{Check --force-single-pool flag}
    F -- True, User wants 1 tube --> G[Run 'Best-Imperfect-Set' v3.2]
    F -- False, Default --> E[Run Auto-N-Pool Logic v4.2]
    E --> H[Start N=1 (Test 1 pool)]
    H --> I{Auto-Healer finds<br>0-clash set?}
    I -- Yes --> J[SUCCESS: Save N=1 files]
    I -- No --> K[Start N=2 (Test 2 pools)]
    K --> L{Auto-Healer finds<br>0-clash set for ALL pools?}
    L -- Yes --> M[SUCCESS: Save N=2 files]
    L -- No --> N[Start N=3...]
    N --> O[...]
    G --> P[Save 'best-available' single pool<br>(may have clashes)]
    J --> Q[End]
    M --> Q
    P --> Q
```

---

## Usage

The script operates in two main modes, which it selects based on the arguments you provide.

### Mode 1: Full Design Pipeline

This is the main, fully-automated mode. It designs, validates, and tails primers from scratch.

**Required Arguments:**
* `--genome`: Path to your reference genome (FASTA format).
* `--gff`: Path to your gene annotation file (GFF format).
* `--target-file`: Path to a simple `.txt` file with one gene name per line.
* `--blast-db`: A prefix for your BLAST database (e.g., "ecoli_db"). The script will automatically build this database from your genome if it doesn't exist.

**Optional Arguments:**
* `--force-single-pool`: (Flag) If provided, it will force the script to find the "best-available" set for a single tube, even if it has clashes. This is only used for `sparse` (non-overlapping) designs.
* `--output-prefix`: (Default: `final_primers`) A prefix for your output files.


#### Example 1: The Default "Auto-Pool" Run (Recommended)
This is the safest, most robust command. The script will auto-detect if your amplicons are sparse or tiled and find the *minimum number of perfect pools* required.

```powershell
python .\design.py --genome "ecoli_genome.fna" --gff "genomic.gff" --target-file "target_genes.txt" --blast-db "ecoli_db" --output-prefix "ecoli_auto_pool_run"
```
* **Output:** Will be `ecoli_auto_pool_run_pool_1.csv`, `..._pool_2.csv`, etc. If a perfect single pool is possible, it will only output one set of files (`..._pool_1.csv`).

#### Example 2: The "Forced Single-Pool" Run
This command forces the script to produce a single-tube solution, *only if* the design is sparse (non-tiled).

```powershell
python .\design.py --genome "ecoli_genome.fna" --gff "genomic.gff" --target-file "target_genes.txt" --blast-db "ecoli_db" --force-single-pool --output-prefix "ecoli_single_pool_run"
```
* **Output:** A single set of files (`.csv`, `.bed`, `.log.txt`). The `.log.txt` will list any clashes it could not resolve. If the design is detected as `tiled`, this flag will be ignored and it will auto-pool anyway.

---

### Mode 2: Tail-Only Utility

This mode bypasses all design and validation. It simply takes existing primer lists and adds your hardcoded sequencing tails.

**Example Command:**
```powershell
python .\design.py --tail-fwd-file "my_fwd_primers.txt" --tail-rev-file "my_rev_primers.txt" --output-prefix "my_tailed_primers"
```

---

## Understanding the Output

The script generates three files for each design or pool:

**1. The `.csv` File (The "Order Sheet")**
This is your main result file. Key columns include:
* **`target_id`**: The gene name for this primer pair.
* **`pair_rank`**: The Primer3 quality score (e.g., `0 (Strategy 0)` is the best primer from the default strategy).
* **`flags`**: **(IMPORTANT)** A direct instruction for lab-bench optimization. A `Low_Tm` or `High_Tm` flag suggests you may need to adjust this primer pair's concentration in your final pool. `OK` means it passed all ideal checks.
* **`fwd_primer_seq` / `rev_primer_seq`**: The 20-mer target-specific sequence.
* **`fwd_primer_tailed` / `rev_primer_tailed`**: The final, 50mer+ sequence (reverse-complemented and tailed), ready to be ordered.
* **`fwd_primer_tm` / `rev_primer_tm`**: The $T_m$ of the *target-specific* sequence.
* **`specificity_hits`**: The BLAST result. `F:1, R:1` is the ideal.

**2. The `.bed` File (The "Visual Check")**
* This is a standard genomics file. Load your reference genome (`.fna`) into a genome browser (like **IGV**), then load this `.bed` file to visually confirm your amplicons are in the correct locations.

**3. The `.log.txt` File (The "Traceability Log")**
* This file contains a complete record of the design process. It lists all targets that failed the initial design and, most importantly, provides the full list of all potential dimer clashes that the "Auto-Healing" logic had to resolve.
* Use this to cross-reference with your sequencing data to troubleshoot failed amplicons.

        </textarea>
        
        <button id="copyButton" class="mt-4 w-full bg-blue-600 text-white font-bold py

