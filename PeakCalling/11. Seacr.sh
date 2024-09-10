## Unnormalized Peak Calling with SEACR ##

# First move IgG control samples to a separate folder
DIR=/path/to/ephemeral/CT_<batch>

mv $DIR/alignment/bedgraph/CT_0*_IgG*_rmDup_bowtie2.fragments.bedgraph $DIR/alignment/bedgraph/histControl

# _seacr.sh is a script to run SEACR
# Run by batch to match to the correct IgG control

#!/bin/bash
#PBS -l select=1:mem=24gb:ncpus=2
#PBS -l walltime=48:00:00
#PBS -N seacr_unnorm_WT_CT_<batch>
#PBS -J 0-5

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate cutnTag

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>

# List your input files
files=$(ls $DIR/alignment/bedgraph/*WT*_rmDup_bowtie2.fragments.bedgraph)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _rmDup_bowtie2.fragments.bedgraph)

# Create path to IgG control with the highest sequencing depth
histControl=$DIR/alignment/bedgraph/histControl/CT_<batch>_WT_IgG_B_rmDup_bowtie2.fragments.bedgraph
cp $histControl .

# Set location of SEACR program
seacr="/path/to/seacr/SEACR_1.3.sh"

#run seacr with built-in normalisation
bash $seacr $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph $histControl norm stringent $DIR/peakCalling/"$sampleID"_rmDup_seacr_control.peaks

#run seacr with built-in normalisation for top 10% peaks
bash $seacr $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph 0.01 norm stringent $DIR/peakCalling/"$sampleID"_rmDup_seacr_top0.01.peaks

#run seacr without normalisation with IgG control
bash $seacr $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph $histControl non stringent $DIR/peakCalling/"$sampleID"_rmDup_seacr_non_control.peaks

#run seacr without normalisation with IgG control in relaxed mode
bash $seacr $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph $histControl non relaxed $DIR/peakCalling/"$sampleID"_rmDup_seacr_non_relaxed.peaks

#run seacr without normalisation without IgG control in stringent mode
bash $seacr $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph 0.01 non stringent $DIR/peakCalling/"$sampleID"_rmDup_seacr_non_top0.01.peaks

#run seacr without normalisation without IgG control in relaxed mode
bash $seacr $DIR/alignment/bedgraph/"$sampleID"_rmDup_bowtie2.fragments.bedgraph 0.01 non relaxed $DIR/peakCalling/"$sampleID"_rmDup_seacr_non_relaxed_top0.01.peaks

