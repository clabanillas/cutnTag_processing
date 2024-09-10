## Assessing replicate reproducibility ##

# To study the reproducibility between replicates and across conditions, the genome is split into 500 bp bins, 
# and a Pearson correlation of the log2-transformed values of read counts in each bin is calculated between replicate datasets.
# Multiple replicates and IgG control datasets can then be displayed in a hierarchically clustered correlation matrix in R

# bin500_bed.sh is a script to split the bed files into 500 bp bins.

#!/bin/bash
#PBS -l select=1:mem=96gb:ncpus=4
#PBS -l walltime=48:00:00
#PBS -N CT_<batch>_binLen500
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate cutnTag

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>/alignment

# List your input files
files=$(ls $DIR/bed/*_rmDup_bowtie2.fragments.bed)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _rmDup_bowtie2.fragments.bed)

# Set the bin length to 500 bp
binLen=500

# Use the midpoint of each fragment to infer which 500 bp bin the fragment belongs to.
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $DIR/bed/"$sampleID"_rmDup_bowtie2.fragments.bed | \
sort -k1,1V -k2,2n | uniq -c | \
awk -v OFS="\t" '{print $2, $3, $1}' | \
sort -k1,1V -k2,2n > $DIR/bed/"$sampleID"_rmDup_bowtie2.fragmentsCount.bin$binLen.bed
