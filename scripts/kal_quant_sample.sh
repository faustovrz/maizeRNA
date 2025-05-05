#!/bin/bash

# Load required modules
conda activate /usr/local/usrapps/maize/sorghum/conda/envs/rnaseq

# Set the path to the kallisto index
index=Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx
# Set the path to the transcripts fasta file
transcripts=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa

# Set the suffix for read 1 and read 2
read1Suffix="_R1_001.fastq.gz"
read2Suffix="_R2_001.fastq.gz"

# Set the path to the fastq files directory
fastqDir=/rsstu/users/r/rrellan/sara/RNAseq/RellanAlvarez

# Get the sample name as a command line argument
sample=$1

# Set the path to the output directory
outputDir=quant_out/${sample}

# Set the paths to the fastq files for the current sample
read1=${fastqDir}/${sample}${read1Suffix}
read2=${fastqDir}/${sample}${read2Suffix}

# Create the output directory
mkdir -p ${outputDir}

# Run kallisto quantification
kallisto quant -i ${index} -o ${outputDir} ${read1} ${read2}
