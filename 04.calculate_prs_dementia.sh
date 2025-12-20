#!/bin/bash
###############################################################################
# Purpose: Calculate Polygenic Risk Scores (PRS) for dementia using PRS-CS
# Description: 
#   This script calculates PRS for two dementia GWAS datasets:
#   1. FinnGen R12 F5 DEMENTIA
#   2. MVP (GCST90475542)
#   
#   This script uses PRS-CS with default settings and automatic phi estimation.
#   It performs the following steps for each dataset:
#   1. Run PRS-CS to estimate SNP effect sizes
#   2. Calculate PRS using PLINK
#
# Software requirements:
#   - PRS-CS (https://github.com/getian107/PRScs)
#   - PLINK 1.9
#
# Input files needed:
#   - 1000 Genomes European LD reference panel
#   - UK Biobank BIM file
#   - Summary statistics files (FinnGen and MVP)
#
# Output files:
#   - PRS effect size files (*.txt)
#   - PLINK profile files (*.profile)
#
###############################################################################

# Set number of threads for parallel processing
N_THREADS=40
export MKL_NUM_THREADS=$N_THREADS      # For Math Kernel Library
export NUMEXPR_NUM_THREADS=$N_THREADS  # For NumExpr
export OMP_NUM_THREADS=$N_THREADS      # For OpenMP

# Set output directory and create if doesn't exist
OUT="/outdir"
mkdir -p "$OUT"

# Paths to reference data and UK Biobank files
ref_1kg_eur="/data/ldblk_1kg_eur"  # 1000 Genomes European LD reference
UKB_bim="/final_filtered"          # UK Biobank BIM file

# Function to process PRS
process_prs() {
    local dataset="$1"
    local SampleSize="$2"
    local Sumstat="$3"
    local SUMDIR="$4"
    
    echo "Processing $dataset"
    
    # Create directories
    mkdir -p "$SUMDIR"
    local OUT_SUBDIR="${OUT}/"
    mkdir -p "$OUT_SUBDIR"
    
    # Run PRScs
    python /PRScs.py \
        --ref_dir "$ref_1kg_eur" \
        --bim_prefix "$UKB_bim" \
        --sst_file "$Sumstat" \
        --n_gwas "$SampleSize" \
        --out_dir "$OUT_SUBDIR" \
        --seed 123
    
    # Combine all chromosome score files
    cat "${OUT_SUBDIR}"*.txt > "${SUMDIR}/score.txt"
    
    # Calculate PRS scores using PLINK
    plink --bfile "$UKB_bim" \
        --score "${SUMDIR}/score.txt" 2 4 6 sum \
        --out "${SUMDIR}/score"
}

# Process finngen_R12_F5_DEMENTIA
process_prs "finngen_R12_F5_DEMENTIA" \
    494845 \
    "/Dementia_Finngen_1.txt" \
    "/sum_res_finngen"

# Process MVP GCST90475542
process_prs "MVP GCST90475542" \
    315668 \
    "/Dementia_MVP_1.txt" \
    "/sum_res_MVP"