#!/bin/bash

# Load required modules
# module load STAR/2.7.7a-foss-2020b

STAR=/usr/local/usrapps/maize/STAR-2.7.10b/bin/Linux_x86_64_static/STAR

# Set the path to the genome reference directory
genomeDir=/rsstu/users/r/rrellan/sara/ref/NAM5_CHR

# Set the path to the annotation file
annotation=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf

# Set the path to the fastq files directory
fastqDir=/rsstu/users/r/rrellan/sara/RNAseq/RellanAlvarez

# Set the path to the output directory
outputdir=counts_out
bamdir=aln_out

# Get the sample name as a command line argument
sample=$1

bamfile=$bamdir/${sample}Aligned.sortedByCoord.out.bam

# Create the output directory
mkdir -p ${outputdir}

# Run STARcount on the current sample's bam files
${STAR} \
     --quantMode GeneCounts \
     --runMode inputAlignmentsFromBAM \
     --inputBAMfile ${bamfile} \
     --runThreadN 8