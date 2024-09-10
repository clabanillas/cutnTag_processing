#!/bin/bash
## Scaling Factor Application ##

# Set directories
DIR=/path/to/ephemeral/CT_<batch>

# Load chromosome sizes file
chromSize=/path/to/genome_indexes/chm13v2.0/chm13v2.0.chrom.sizes

# Activate the environment
module load anaconda3/personal
source activate cutnTag

# Load processing information (including scaling factor) into an array
arr_csv=()
while IFS= read -r line
do
  arr_csv+=("$line")
done < $DIR/kmet/Proc_Info_CT_<batch>.csv

# Apply scaling factor and generate normalized bedgraph files
COUNTER=0
for i in $DIR/alignment/bed/*.fragments.bed
do 
  COUNTER=$((COUNTER+1))
  myline=${arr_csv[$COUNTER]}
  IFS=',' read -r -a myline2 <<< "$myline"
  scale_factor=${myline2[-1]}  # Extract scale factor
  sampleID=$(basename "$i" .fragments.bed)

  # Generate normalized bedgraph using bedtools
  bedtools genomecov -bg -scale $scale_factor -i $i -g $chromSize \
  > $DIR/alignment/bedgraph/normalised/${sampleID}.normalised.bedgraph
done
