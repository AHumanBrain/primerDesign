Multiplex PCR Primer Designer

A command-line tool for designing specific and compatible primers for multiplex PCR, with a focus on preparing amplicons for next-generation sequencing.

Overview

This tool automates the process of designing primers for targeted sequencing. It takes a reference genome and a set of target regions (genes, coordinates, etc.) and uses a combination of Primer3 and BLAST to find optimal primer pairs. The core logic is built to handle specificity checks and multiplex compatibility, ensuring the designed primers work well together in a single reaction.

Key Features

The development is planned in phases:

Core Primer Design:

Designs standard primer pairs for single genomic targets using the robust Primer3 engine.

Accepts genomes (FASTA) and target regions (gene names via GFF, or coordinates via BED) as input.

Outputs results in user-friendly formats, including a .csv for primer details and a .bed file for visualizing amplicons in a genome browser.

Specificity & Multiplexing:

Specificity Guarantee: Integrates with local BLAST to ensure designed primers are specific to the target region and will not produce off-target amplicons.

Multiplex Compatibility: Performs in silico checks for primer-dimer formation (self-dimers and cross-dimers) across the entire primer pool.

Automated Pooling: Groups primers into compatible pools to maximize multiplexing efficiency.

Usability & Advanced Features:

Adapter Tailing: Optionally appends sequencing adapters (e.g., for Illumina, PacBio) to the 5' end of designed primers.

Configuration Files: Allows for easy management of primer design parameters (Tm, size, GC content, etc.) through a simple configuration file.

Detailed Logging: Provides clear feedback on the design process, including reasons for primer rejection.

Installation

Clone the repository:

git clone [https://github.com/your-username/multiplex-primer-designer.git](https://github.com/your-username/multiplex-primer-designer.git)
cd multiplex-primer-designer


Create and activate a virtual environment (recommended):

python3 -m venv venv
source venv/bin/activate


Install dependencies:

This tool requires a local installation of NCBI BLAST+. You can find installation instructions here.

Install the required Python packages:

pip install -r requirements.txt


Usage

The basic command structure will be as follows. First, create a BLAST database for your reference genome:

makeblastdb -in path/to/your/ecoli_genome.fasta -dbtype nucl -out ecoli_db


Then, run the primer design script:

python design_primers.py \
    --genome path/to/your/ecoli_genome.fasta \
    --gff path/to/your/ecoli_genes.gff \
    --target-id gyrA \
    --blast-db ecoli_db \
    --output-prefix ecoli_gyrA_primers


This will generate ecoli_gyrA_primers.csv and ecoli_gyrA_primers.bed.

This project is currently under active development.
