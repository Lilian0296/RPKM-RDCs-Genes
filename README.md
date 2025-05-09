# RPKM-RDCs-Genes

Calculation of Transcribed Regions Within RDCs from GRO-seq data

# Overview
This script is designed to calculate reads per kilobase per millions (RPKM) in RDC regions.

# Mode 

You can assign a mode to different calculation

- **Gene RPKM:**  Calculates RPKM for all genes (all reads)
- **RDC RPKM:** Calculates RPKM for genes overlapping with RDCs (considering gene strandness). For RDCs without any overlapping genes, the calculation is based on the number of reads within the regions, without considering strandness.

# Set up the environment

```bash
# bash
module load R/4.4.3-GCCcore-14.1.0
module load SAMtools/1.20-GCC-14.1.0 
module load Subread/2.0.6-GCC-14.1.0

```

```R
# R
install.packages("readr")
install.packages("dplyr")
install.packages("optparse")

```

# Running the Script
## Input Parameters:
```bash
-o / --work_dir: Working directory. Specifies the main working directory where output files will be saved or accessed.
-r / --genome: Genome reference file (Such as hg38_refGene.bed)
-b / --bam_files: BAM files folder. Path to the directory containing the input BAM files used for read mapping and quantification.
-i / --rdc_input: RDC BED files folder (type: character). Path to the directory containing the input BED files for RDC regions.
```

# Final output

- Gene RPKM: CSV files with RPKM values.
- RDC RPKM: For RDCs that overlap with genes, the file names start with “Gene_“. For RDCs that do not overlap with genes, the file names start with “RDC_“.



# Example

```bash
Rscript 0_RDC_rpkm_V1.R -o ~/u2os_res/ -r hg38_refGene.bed -b ~/gro-seq-u2os-hg38/alignment/ -i ~/u2os/ 

```

