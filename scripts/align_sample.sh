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
# this only writes Log.out 
outputDir=aln_out

# Set the suffix for read 1 and read 2
read1Suffix="_R1_001.fastq.gz"
read2Suffix="_R2_001.fastq.gz"

# Get the sample name as a command line argument
sample=$1

# Set the paths to the fastq files for the current sample
read1=${fastqDir}/${sample}${read1Suffix}
read2=${fastqDir}/${sample}${read2Suffix}

# Set the output file name for the current sample
outputFile=${outputDir}/${sample}

# Create the output directory
mkdir -p ${outputDir}

# Run STAR alignment on the current sample's fastq files and indexed genome
${STAR} \
     --runThreadN 8 \
     --genomeDir ${genomeDir} \
     --sjdbGTFfile ${annotation} \
     --readFilesCommand zcat \
     --readFilesIn ${read1} ${read2} \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${outputFile} \
     --outReadsUnmapped Fastx