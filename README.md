# Multiplex PCR Primer Designer CLI

This is a command-line tool for designing specific, adapter-tailed primer pairs for multiplex PCR. The tool is built in Python and leverages `primer3-py` for primer design and a local NCBI BLAST+ installation for specificity checking.

It features two main modes of operation:

1.  **Full Design Mode:** An end-to-end pipeline that takes a genome and a list of target genes, designs specific primers, and adds sequencing adapter tails.
2.  **Tail-Only Mode:** A utility to quickly add adapter tails to pre-existing lists of forward and reverse primers.

## Features

* **Primer3 Integration:** Designs thermodynamically optimal primers based on your sequence.
* **Automated Specificity Check:** Automatically runs each candidate primer pair against a local BLAST database and filters for pairs where both primers have exactly one perfect hit in the genome.
* **Robust File Handling:** Automatically builds the required BLAST database from your genome, cleaning FASTA headers to prevent common errors.
* **Custom Tailing Logic:** Hardcoded with specific adapter sequences, which are added to the **reverse complement** of the designed primers.
* **Dual Output:**
    * `.csv`: A detailed report with all specific primer pairs, their properties (Tm, size), and the final tailed sequences.
    * `.bed`: A browser-ready file to visualize the location of your amplicons in a genome browser like IGV.

## Setup

### 1. Clone the Repository

```bash
git clone [https://github.com/AHumanBrain/primer_designer.git](https://github.com/AHumanBrain/primer_designer.git)
cd primer_designer
```

### 2. Install NCBI BLAST+

This tool **requires** a local installation of NCBI BLAST+.

1.  **Download:** Get the installer from the [NCBI BLAST+ Download Page](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/). Choose the `x64-win64.exe` file.
2.  **Install:** Run the installer. **Crucially, ensure you check the box that says "Add BLAST to the system PATH".**
3.  **Verify:** Close and re-open your terminal/PowerShell. Type `blastn -version`. You should see a version number printed.

### 3. Install Python Dependencies

Make sure you are in the project's root directory (where `requirements.txt` is located).

```powershell
pip install -r requirements.txt
```

## Usage

The script operates in one of two modes, determined by the arguments you provide.

### Mode 1: Full Design Pipeline

This mode designs primers from scratch. It requires a genome, a gene annotation file, a target list, and a name for the BLAST database.

**Note:** The first time you run this command on a new genome, the script will automatically create the necessary BLAST database files (e..g, `ecoli_db.nin`, `ecoli_db.nhr`, etc.). This may take a few moments. Subsequent runs will be much faster as they will re-use the existing database.

#### Example Command:

```powershell
python .\design.py --genome "path\to\ecoli_genome.fna" ^
                  --gff "path\to\ecoli_annotations.gff" ^
                  --target-file "path\to\my_targets.txt" ^
                  --blast-db "ecoli_db" ^
                  --output-prefix "ecoli_run_1"
```

* `--genome`: Path to your genome FASTA file.
* `--gff`: Path to your gene annotation file (GFF/GFF3).
* `--target-file`: A simple `.txt` file with one gene name (matching the GFF) per line.
* `--blast-db`: The desired prefix for your BLAST database (e.g., "ecoli\_db").
* `--output-prefix`: The base name for your output files (e.g., `ecoli_run_1.csv` and `ecoli_run_1.bed`).

### Mode 2: Tail-Only Utility

This mode skips all design and BLAST steps. It simply reads your primer lists, applies the hardcoded tailing logic, and saves the output.

First, create your primer files (e.g., `my_fwd_primers.txt` and `my_rev_primers.txt`) with one primer sequence per line.

#### Example Command:

```powershell
python .\design.py --tail-fwd-file "my_fwd_primers.txt" ^
                  --tail-rev-file "my_rev_primers.txt" ^
                  --output-prefix "my_tailed_primers"
```

* `--tail-fwd-file`: Path to your file of forward primers.
* `--tail-rev-file`: Path to your file of reverse primers.
* `--output-prefix`: The base name for your output `.csv` file.

## File Descriptions

* `design.py`: The main Python script.
* `requirements.txt`: A list of required Python libraries.
* `.gitignore`: Standard Python gitignore file.
* `README.md`: This file.
        
