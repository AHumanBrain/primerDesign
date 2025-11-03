Grant-Ready Protocol Copier
The text box below contains the raw Markdown for the protocol. Use the button to copy it.

# Protocol: High-Plex Amplicon Library Preparation for NGS

## I. Summary
This protocol details a two-stage, high-fidelity PCR method for the targeted amplification of N-plex primer panels (designed *in silico*) for next-generation sequencing on Illumina platforms. The workflow consists of: 1) A limited-cycle target enrichment PCR, 2) An enzymatic and size-selection cleanup, and 3) A universal library amplification PCR.

## II. Reagents & Equipment
* **DNA:** 10-100 ng high-quality genomic DNA (gDNA) per sample.
* **Enzymes:**
    * High-fidelity, Hot-Start DNA Polymerase (e.g., Q5, KAPA HiFi, Phusion).
    * Exonuclease I (e.g., NEB).
* **Reagents:**
    * 5X Polymerase Reaction Buffer.
    * 10 mM dNTP mix.
    * Nuclease-free water.
    * (Optional) 5M Betaine (for GC-rich targets).
* **Primers (from *in silico* design):**
    * **Custom Primer Pool(s):** Lyophilized, desalted primers from *in silico* design, reconstituted to 100 µM in 1X TE Buffer. A working stock (e.g., 2 µM combined-pair concentration) should be prepared.
    * **Universal Primers:** `FWD_TAIL` and `REV_TAIL` sequences, ordered as standard oligos.
* **Cleanup:**
    * SPRI Beads (e.g., AMPure XP or equivalent).
    * Freshly prepared 80% Ethanol.
    * 1X TE Buffer (or low-EDTA TE).
* **Equipment:**
    * Thermocycler.
    * Magnetic stand for 96-well plates or tubes.
    * Fluorometric quantitation system (e.g., Qubit).
    * Capillary electrophoresis system (e.g., Agilent Bioanalyzer or TapeStation).

---

## III. Library Preparation Workflow

### A. Stage 1: Target Enrichment PCR (Multiplex)

This step enriches for the specific target regions from the gDNA. If using a multi-pool design (e.g., `_pool_A`, `_pool_B`), a separate PCR must be run for each pool.

**1. Reaction Recipe (25 µL total volume)**
| Component | Stock Conc. | Final Conc. | Volume (µL) |
| :--- | :---: | :---: | :---: |
| 5X HiFi Buffer | 5X | 1X | 5.0 |
| 10 mM dNTPs | 10 mM | 200 µM | 0.5 |
| Custom Primer Pool (e.g., Pool A) | 2 µM | **50 nM** | 0.625 |
| gDNA | 10 ng/µL | 10-100 ng | 1.0 - 10.0 |
| (Optional) Betaine | 5 M | 1 M | 5.0 |
| Hot-Start HiFi Polymerase | 2 U/µL | 0.02 U/µL | 0.25 |
| Nuclease-free Water | - | - | Up to 25.0 |

**2. Thermalcycling Conditions (Stage 1)**
| Step | Temperature | Time | Cycles |
| :--- | :---: | :---: | :---: |
| 1. Initial Denaturation | 98°C | 3 min | 1x |
| 2. Denaturation | 98°C | 20 sec | \multirow{2}{*}{**10-15**} |
| 3. **Anneal / Extend** | **65°C** | **3 min** | |
| 4. Final Extension | 72°C | 5 min | 1x |
| 5. Hold | 4°C | $\infty$ | |

**3. (If Multi-Pool) Combine Products:**
* Following Stage 1, combine the full 25 µL reactions from Pool A and Pool B (and C, D...) for each sample into a single tube. Proceed to cleanup with this combined 50-100 µL volume.

### B. Stage 2: Enzymatic & Size-Selection Cleanup

This crucial step removes unused primers (ssDNA) and primer-dimers (dsDNA) from the reaction.

1.  **Exonuclease I Digestion (Removes ssDNA Primers):**
    * Add 2 µL of Exonuclease I (NEB, 20 U/µL) directly to the 50 µL combined PCR product.
    * Incubate at **37°C for 30 minutes**.
    * Incubate at **80°C for 20 minutes** to heat-inactivate the exonuclease.

2.  **SPRI Bead Cleanup (Removes Primer-Dimers):**
    * Warm SPRI beads to room temperature (RT) for 15 min.
    * Add **0.8X** volume of SPRI beads to the reaction (e.g., 40 µL of beads for a 50 µL reaction). *This ratio is a critical optimization point.*
    * Pipette mix 10 times and incubate at RT for 5 min to bind DNA.
    * Place tube on magnetic stand for 2 min until clear.
    * Aspirate and discard the supernatant (which contains primer-dimers).
    * Wash the beads twice with 200 µL of 80% Ethanol. Do *not* resuspend.
    * Air dry beads on the magnet for 3-5 min.
    * Remove from magnet and elute DNA in 22 µL of 1X TE Buffer.

### C. Stage 3: Universal Library Amplification PCR

This step adds the full sequencing adapters and barcodes (if applicable) using the tails as priming sites.

**1. Reaction Recipe (50 µL total volume)**
| Component | Stock Conc. | Final Conc. | Volume (µL) |
| :--- | :---: | :---: | :---: |
| 5X HiFi Buffer | 5X | 1X | 10.0 |
| 10 mM dNTPs | 10 mM | 200 µM | 1.0 |
| Universal FWD Primer (FWD_TAIL) | 10 µM | 0.5 µM | 2.5 |
| Universal REV Primer (REV_TAIL) | 10 µM | 0.5 µM | 2.5 |
| Cleaned DNA (from Step B) | - | - | 20.0 |
| Hot-Start HiFi Polymerase | 2 U/µL | 0.02 U/µL | 0.5 |
| Nuclease-free Water | - | - | 13.5 |

**2. Thermalcycling Conditions (Stage 2)**
| Step | Temperature | Time | Cycles |
| :--- | :---: | :---: | :---: |
| 1. Initial Denaturation | 98°C | 3 min | 1x |
| 2. Denaturation | 98°C | 20 sec | \multirow{2}{*}{**15-20**} |
| 3. Anneal / Extend | 70°C | 1 min | |
| 4. Final Extension | 72°C | 5 min | 1x |
| 5. Hold | 4°C | $\infty$ | |

*Note: A second 0.8X SPRI cleanup is typically performed after this step to purify the final library before quantification.*

---

## IV. Library QC and Normalization
1.  **Quantification:** Quantify the final library using a fluorometric method (e.g., Qubit dsDNA HS Assay).
2.  **Sizing:** Assess the library size distribution and confirm the absence of primer-dimers using a capillary electrophoresis system (e.g., Agilent TapeStation or Bioanalyzer). The library should appear as a clean peak in the 200-300bp range (150-250bp amplicons + adapter tails).

---

## V. Key Optimization Parameters
This protocol is a starting point. For novel panels, optimization is required to ensure balanced amplification.
* **Stage 1 Primer Concentration:** The 50 nM concentration is a "starvation" level. If yield is low, this can be increased (e.g., to 75-100 nM) at the risk of higher dimer formation.
* **Stage 1 Anneal/Extend (A/E) Temp:** The 65°C A/E temperature is stringent. If amplicons with low $T_m$ (per the design log) are under-represented, this temperature can be lowered (e.g., to 60-62°C).
* **Stage 1 Cycle Count:** This is the most effective variable for increasing yield. For low-input gDNA, increase from 15 to 18 or 20 cycles.
* **Betaine:** For high-GC panels (per the design log), add 1M Betaine to the Stage 1 PCR to improve the amplification of GC-rich, "knotted" targets.
* **SPRI Bead Ratio:** The 0.8X ratio is designed to cut off fragments <100bp. This ratio may need to be adjusted based on the observed size of primer-dimers.
* **Primer Pool Balancing:** For final optimization, use sequencing coverage data to "re-balance" the primer pool. Increase the concentration of primers for low-coverage amplicons and decrease the concentration for high-coverage amplicons.

        
 Copy Raw Markdown to Clipboard
