# This series of scripts can be used to build a Transcription Start Site heatmap of samples .bam files over a genome .gtf


## 1.Generate raw bigwig 


# bam_to_bigwig is a script to convert .bam files to .bw (BigWig) 

#!/bin/bash
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N bam_to_bigwig_CT_<batch>
#PBS -J 0-11

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
# Load python environment with deeptools in it
source activate cutnTag2

# Set directories
DIR=/path/to/ephemeral/CT_<batch>
OUTDIR=/path/to/ephemeral/CT_<batch>/heatmaps

# List your input files
files=$(ls $DIR/alignment/bam/*_bowtie2.rmDup.mapped.bam)
arr=($files)
sample=${arr[$tid]}

# Create sampleID, common for matched _1 and _2 files
sampleID=$(basename $sample _bowtie2.rmDup.mapped.bam)

# Sort by order name
samtools sort -o $OUTDIR/"$sampleID"_bowtie2.rmDup.osorted.mapped.bam $DIR/alignment/bam/"$sampleID"_bowtie2.rmDup.mapped.bam

# Index the sorted bam file
samtools index $OUTDIR/"$sampleID"_bowtie2.rmDup.osorted.mapped.bam

# Generate bigWig from BAM file
bamCoverage -b $OUTDIR/"$sampleID"_bowtie2.rmDup.osorted.mapped.bam -o $OUTDIR/"$sampleID"_raw.bw





## 2. Generate compute matrix


## TSS_heatmap.sh is a script used to generate a heatmap over transcription start sites

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N TSS_heatmap_CT_<batch>

module load anaconda3/personal
source activate cutnTag2

DIR=/path/to/ephemeral/CT_<batch>/alignment/bigwig
genes=/path/to/gene/directory/genes
OUTDIR=/path/to/ephemeral/CT_<batch>/heatmaps/TSS

cores=8

computeMatrix scale-regions -S $DIR/CT_<batch>_<sample1>_raw.bw \
                               $DIR/CT_<batch>_<sample2>_raw.bw \
                              -R $genes/T2T-CHM13v2.0_genomic_modified.gtf \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $OUTDIR/CT_<batch>_matrix_T2T_TSS.mat.gz -p $cores

plotHeatmap -m $OUTDIR/CT_<batch>_matrix_T2T_TSS.mat.gz -out $OUTDIR/CT_<batch>_T2T_TSS_heatmap.png --sortUsing sum









