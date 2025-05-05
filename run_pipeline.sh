#!/bin/bash

# Maize RNA-seq Analysis Pipeline
# This script orchestrates the full analysis workflow

echo "Starting RNA-seq analysis pipeline..."

# 1. Index the genome with STAR
echo "Step 1: Indexing genome with STAR..."
bsub < batch/q_STAR_index_genome_chr.sh

# Wait for the indexing job to complete
# In a real environment, you would check job status or use job dependencies
sleep 2

# 2. Index the transcriptome with Kallisto
echo "Step 2: Indexing transcriptome with Kallisto..."
bsub < batch/q_kallisto_idx.sh

# 3. Submit alignment jobs for all samples
echo "Step 3: Aligning reads with STAR..."
bsub < batch/q_STAR_align.sh sample_list.txt

# 4. Quantify expression with Kallisto
echo "Step 4: Quantifying expression with Kallisto..."
bsub < batch/q_kal_quant_sample.sh sample_list.txt

# 5. Run statistical analysis
echo "Step 5: Performing expression analysis..."
# This would typically be run once all previous jobs have completed
Rscript scripts/kallisto_stats.R

echo "Pipeline submitted. Check job status with bjobs."