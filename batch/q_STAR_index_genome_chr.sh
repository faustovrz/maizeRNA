#!/bin/bash
#BSUB -J star_index
#BSUB -q sara
#BSUB -n 24
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=128GB]"
#BSUB -W 24:00
#BSUB -o %J.out
#BSUB -e %J.err

# Execute the STAR indexing script
../scripts/STAR_index_genome_chr.sh