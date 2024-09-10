# This script takes a peak file and calculates read counts over peaks vs over rest of genome (readout of signal-to-noise ratio)

#!/bin/bash
#PBS -l select=1:mem=128gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N peak_heatmap_CT_<batch>
#PBS -J 0-11

tid=$PBS_ARRAY_INDEX

module load anaconda3/personal
source activate cutnTag2

DIR=/path/to/ephemeral/CT_<batch>/alignment/bigwig
PEAKDIR=/path/to/ephemeral/CT_<batch>/peakCalling
OUTDIR=/path/to/ephemeral/CT_<batch>/heatmaps/peak

files=$(ls $DIR/*_raw.bw)
arr=($files)
sample=${arr[$tid]}

sampleID=$(echo `basename $sample _raw.bw`)

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $PEAKDIR/"$sampleID"_seacr_control.peaks.stringent.bed >$PEAKDIR/"$sampleID"_seacr_control.peaks.summitRegion.bed

cores=8
computeMatrix reference-point -S $DIR/"$sampleID"_raw.bw \
              -R $PEAKDIR/"$sampleID"_seacr_control.peaks.summitRegion.bed \
              --skipZeros -o $OUTDIR/"$sampleID"_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $OUTDIR/"$sampleID"_SEACR.mat.gz -out $OUTDIR/"$sampleID"_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "$sampleID"
