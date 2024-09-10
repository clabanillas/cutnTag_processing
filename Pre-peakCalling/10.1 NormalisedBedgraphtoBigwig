## Normalized BigWig Conversion ##

#!/bin/bash
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N bigwig_norm_<batch>
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate ucsc-tools

# Set directories
DIR=/path/to/ephemeral/CT_<batch>/alignment/bedgraph/normalised
OUTDIR=/path/to/ephemeral/CT_<batch>/peakCalling/normalised/bigwig
chromSizes=/path/to/genome_indexes/chm13v2.0/chm13v2.0.chrom.sizes

# List your input files
files=$(ls $DIR/*_rmDup_bowtie2.normalised.bedgraph)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _rmDup_bowtie2.normalised.bedgraph)

# Convert bedgraph to BigWig
bedGraphToBigWig $DIR/"$sampleID"_rmDup_bowtie2.normalised.bedgraph $chromSizes $OUTDIR/"$sampleID"_normalised_rmDup.bw
