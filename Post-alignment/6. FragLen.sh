## Assess mapped fragment size distribution ##

# FragLen.sh is a script that runs samtools on a list of .sam samples to assess mapped fragment size distribution

#!/bin/bash
#PBS -l select=5:mem=124gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N CT_<batch>_fraglength
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate cutnTag

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>/alignment

# List your input files
files=$(ls $DIR/sam/*_bowtie2.sam)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _bowtie2.sam)

# Extract the 9th column from the alignment sam file, which is the fragment length
samtools view -F 0x04 $DIR/sam/"$sampleID"_bowtie2.sam | \
awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > $DIR/sam/fragmentLen/"$sampleID"_fragmentLen.txt
