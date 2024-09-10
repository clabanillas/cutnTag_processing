## removing duplicates ##

#_remdup_picard.sh is a script that runs picard tools on a list of .sam samples to (1)sort (2)mark and (3)remove duplicates

#!/bin/bash
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N remdup_picard_<batch>
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate cutnTag

# Set working directories
DIR=/path/to/CT_<batch>/alignment
PICARD=/path/to/picard/picard-2.27.4-hdfd78af_0/share/picard-2.27.4-0/picard.jar
ephem=/path/to/ephemeral

# List your input files
files=$(ls $DIR/sam/*_bowtie2.sam)
arr=($files)
sample=${arr[$tid]}

# Create sampleID
sampleID=$(basename $sample _bowtie2.sam)

# Sort by coordinate using Picard's SortSam tool
java -jar $PICARD SortSam I=$DIR/sam/"$sampleID"_bowtie2.sam O=$DIR/sam/"$sampleID"_bowtie2.sorted.sam SORT_ORDER=coordinate

# Mark duplicates using Picard's MarkDuplicates tool
java -jar $PICARD MarkDuplicates I=$DIR/sam/"$sampleID"_bowtie2.sorted.sam O=$DIR/removeDuplicate/"$sampleID"_bowtie2.sorted.dupMarked.sam METRICS_FILE=$DIR/removeDuplicate/picard_summary/"$sampleID"_picard.dupMark.txt TMP_DIR=$ephem/TMP

# Remove duplicates
java -jar $PICARD MarkDuplicates I=$DIR/sam/"$sampleID"_bowtie2.sorted.sam O=$DIR/removeDuplicate/"$sampleID"_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$DIR/removeDuplicate/picard_summary/"$sampleID"_picard.rmDup.txt TMP_DIR=$ephem/TMP
