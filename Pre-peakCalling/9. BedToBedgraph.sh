## Generating unnormalized bedgraphs ##

# BedToBedgraph.sh is a script to convert bed files to bedgraph format

#!/bin/bash
#PBS -l select=1:mem=94gb:ncpus=4
#PBS -l walltime=48:00:00
#PBS -N CT_<batch>_bed_to_bedgraph
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate cutnTag

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>

# List your input files
files=$(ls $DIR/alignment/bed/*_rmDup_bowtie2.fragments.bed)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _rmDup_bowtie2.fragments.bed)

# Path to chromosome sizes file
chromSizes=/path/to/genome_indexes/chm13v2.0/chm13v2.0.chrom.sizes
cp $chromSizes .

# Generate bedgraph from fragments.bed files
bedtools genomecov -bg -i $sample -g chm13v2.0.chrom.sizes > $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph
