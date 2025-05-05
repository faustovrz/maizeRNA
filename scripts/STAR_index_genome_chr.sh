#!/bin/bash

# Load required modules
#module load STAR/2.7.7a-foss-2020b

STAR=/usr/local/usrapps/maize/STAR-2.7.10b/bin/Linux_x86_64_static/STAR

# Set the path to the genome reference directory
# in this directory the index SAindex will be written
genomeDir=/rsstu/users/r/rrellan/sara/ref/NAM5_CHR

# Set the path to the annotation file
annotation=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.56.gtf

# Set the path to the output directory
# this only writes Log.out 
outputDir=aln_out

# Set the number of threads
threads=24

# Create the output directory
mkdir -p ${outputDir}

# Create the STAR index
${STAR} \
  --runThreadN ${threads} \
  --runMode genomeGenerate \
  --genomeDir ${genomeDir} \
  --genomeFastaFiles ${genomeDir}/*.fa \
  --sjdbGTFfile ${annotation} \
  --sjdbOverhang 100 \
  --outFileNamePrefix ${outputDir}/