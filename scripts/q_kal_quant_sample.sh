#!/bin/bash

# Set the path to the sample list file
sampleList=$1
index=Zm-B73-REFERENCE-NAM-5.0.cdna.all.idx
transcripts=/rsstu/users/r/rrellan/sara/ref/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cdna.all.fa

# Loop over the sample names in the sample list file
while read sample; do
  # Construct the bsub command to run the alignment job for the current sample
  # par="-q sara -n 12 -R 'span[hosts=1]' -R 'rusage[mem=32GB]' -W 4:00 -o %J.stdout -e %J.stderr" 
  par="-n 1 -W 4:00 -o %J.stdout -e %J.stderr" 
  echo $sample 
  echo $par
  # Submit the bsub command to the job scheduler
  bsub -R 'rusage[mem=2GB]' ${par} ../scripts/kal_quant_sample.sh ${sample}
done < ${sampleList}
