## Generating unnormalized bigWigs ##

# _bedgraph_to_bigwig.sh is a script to convert .bedgraph files to .bw (bigWig) format

#!/bin/bash
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N unnorm_bigwig_<batch>
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate ucsc-tools

# Set directories
DIR=/path/to/ephemeral/CT_<batch>/alignment/bedgraph
OUTDIR=/path/to/ephemeral/CT_<batch>/peakCalling/bigwig
chromSizes=/path/to/genome_indexes/chm13v2.0/chm13v2.0.chrom.sizes

# List your input files
files=$(ls $DIR/*_rmDup_bowtie2.fragments.bedgraph)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _rmDup_bowtie2.fragments.bedgraph)

# Convert bedgraph to bigWig
bedGraphToBigWig $DIR/"$sampleID"_rmDup_bowtie2.fragments.bedgraph $chromSizes $OUTDIR/"$sampleID"_rmDup.bw
