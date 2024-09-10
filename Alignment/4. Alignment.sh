## alignment ##


# alignment.sh is a script that indexes and runs Bowtie2 v2.5.1 alignment on a list of samples for read 1 and read 2


#!/bin/bash
#PBS -l select=1:mem=64gb:ncpus=10
#PBS -l walltime=48:00:00
#PBS -N alignment_<batch>
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Load required module
module load anaconda3/personal
source activate bowtie2

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>/fastp
GENOMEDIR=/path/to/project/genome_indexes/chm13v2.0
TEMPDIR=/path/to/ephemeral/CT_<batch>

# List your input files
files=$(ls $DIR/*_1.trimmed.fastq.gz)
arr=($files)
sample=${arr[$tid]}

# Create sampleID, common for matched _1 and _2 files
sampleID=$(basename $sample _1.trimmed.fastq.gz)

# Create read 1 and read 2 input variables and samplename for output
in1=$DIR/"$sampleID"_1.trimmed.fastq.gz
in2=$DIR/"$sampleID"_2.trimmed.fastq.gz

# Run alignment with Bowtie2 v2.5.1
cores=10

# PE + trimmed option
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p $cores -x $GENOMEDIR/chm13v2.0 -1 $in1 -2 $in2 -S $TEMPDIR/alignment/sam/"$sampleID"_bowtie2.sam &> $TEMPDIR/alignment/sam/bowtie2_summary/"$sampleID"_bowtie2.txt


## TERMINAL ##

# Check if files were created
ls -lh $DIR/alignment/sam/*.sam 

# Generate quick alignment report
cd $DIR/alignment/sam/bowtie2_summary
module load anaconda3/personal
source activate cutnTag
multiqc .
