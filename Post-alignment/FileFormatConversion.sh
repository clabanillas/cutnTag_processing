#!/bin/bash
## File format conversion from rmDup files ##

# fileFormatConversion.sh is a script adapted from 131022_rmDup_bam.sh and 080922_bam_bed_samtools.sh to run file conversion

#PBS -l select=1:mem=100gb:ncpus=8 
#PBS -l walltime=48:00:00
#PBS -N CT_<batch>_fileFormatConversion
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate cutnTag

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>/alignment

# List your input files
files=$(ls $DIR/removeDuplicate/*_bowtie2.sorted.rmDup.sam)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _bowtie2.sorted.rmDup.sam)

# Step 1: Filter and keep the mapped read pairs in .bam format
samtools view -bS -F 0x04 $DIR/removeDuplicate/"$sampleID"_bowtie2.sorted.rmDup.sam > $DIR/bam/"$sampleID"_bowtie2.sorted.rmDup.mapped.bam

# Step 2: Sort by read name if rmDup files are used
# NOTE: .bam files need to be sorted by read name (-n) to convert to .bed file
samtools sort -o $DIR/bam/"$sampleID"_bowtie2.rmDup.mapped.bam -n $DIR/bam/"$sampleID"_bowtie2.sorted.rmDup.mapped.bam

# Step 3: Convert .bam into .bed file format
bedtools bamtobed -i $DIR/bam/"$sampleID"_bowtie2.rmDup.mapped.bam -bedpe > $DIR/bed/"$sampleID"_rmDup_bowtie2.bed

# Step 4: Keep read pairs that are on the same chromosome and fragment length less than 1000bp
awk '$1==$4 && $6-$2 < 1000 {print $0}' $DIR/bed/"$sampleID"_rmDup_bowtie2.bed > $DIR/bed/"$sampleID"_rmDup_bowtie2.clean.bed

# Step 5: Extract the fragment-related columns
cut -f 1,2,6 $DIR/bed/"$sampleID"_rmDup_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n > $DIR/bed/"$sampleID"_rmDup_bowtie2.fragments.bed
