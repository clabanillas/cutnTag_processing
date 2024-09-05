```bash

## fastqc ##


# fastqc.sh is a script that runs fastqc analysis on a list of samples within a folder
# Remember to customize PBS -N and PBS -J 0-N 

#!/bin/bash
#PBS -l select=1:mem=15gb:ncpus=1
#PBS -l walltime=02:00:00
#PBS -N fastqc_<batch>
#PBS -J 0-31

tid=$PBS_ARRAY_INDEX

# Set working directory
DIR=/path/to/ephemeral/CT_<batch>

# Load required module
module load anaconda3/personal
source activate cutnTag

# List your input files
files=$(ls $DIR/*/*/*/*/*.fq.gz)
arr=($files)
sample=${arr[$tid]}

# Run FastQC command
fastqc $sample -d . -o .

# Move all output to output folder
mv * $DIR/fastqc


## multiqc ##


# multiqc.sh is a script that combines FastQC analysis on a list of samples within a folder
# Remember to customize PBS -N and PBS -J 0-N 

#!/bin/bash
#PBS -l select=1:mem=15gb:ncpus=1
#PBS -l walltime=02:00:00
#PBS -N multiqc_<batch>

# Set working directory
DIR=/path/to/ephemeral/CT_<batch>/fastqc 

# Load required module
module load anaconda3/personal
source activate cutnTag

# Run MultiQC command
multiqc $DIR

# Move all output to output folder
mv * $DIR/multiqc


## adapter trimming ##

# fastp.sh uses the fastp program to cut adapters from sequencing data

#!/bin/bash
#PBS -l select=1:mem=30gb:ncpus=6
#PBS -l walltime=02:00:00
#PBS -N fastp_<batch>
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Set working directory
DIR=/path/to/ephemeral/CT_<batch>
# Set output directory
OUTDIR=$DIR/fastp

# Load required modules
module load anaconda3/personal
source activate cutnTag

# List your input files
files=$(ls $DIR/*/*/*/*/*_1.fq.gz)
arr=($files)
sample=${arr[$tid]}

# Create sampleID, common for matched _1 and _2 files 
sampleID=$(basename $sample _1.fq.gz)

# Create working directory where the files are located
WORKDIR=$(dirname $sample)

# Running fastp on paired-end reads:
fastp -i $DIR/*/*/*/*/"$sampleID"_1.fq.gz -I $DIR/*/*/*/*/"$sampleID"_2.fq.gz -o "$sampleID"_1.trimmed.fastq.gz -O "$sampleID"_2.trimmed.fastq.gz --detect_adapter_for_pe -l 20 -j "$sampleID".fastp.json -h "$sampleID".fastp.html

# Move all output files to output folder
mv * $OUTDIR


## multiqc of trimming reports ##

#!/bin/bash
#PBS -l select=1:mem=15gb:ncpus=1
#PBS -l walltime=02:00:00
#PBS -N multiqc_fastp_<batch>

# Set working directory
DIR=/path/to/ephemeral/CT_<batch>/fastp

# Load required module
module load anaconda3/personal
source activate cutnTag

# Run MultiQC command
multiqc $DIR/*

# Move all output to output folder
mv * $DIR/multiqc


## fastqc of trimmed reads ##


# fastqc_trimmed.sh is a script that runs fastqc analysis on a list of samples within a folder
# Remember to customize PBS -N and PBS -J 0-N 

#!/bin/bash
#PBS -l select=1:mem=15gb:ncpus=1
#PBS -l walltime=02:00:00
#PBS -N fastqc_trimmed_<batch>
#PBS -J 0-31

tid=$PBS_ARRAY_INDEX

# Set working directory
DIR=/path/to/ephemeral/CT_<batch>

# Load required module
module load anaconda3/personal
source activate cutnTag

# List your input files
files=$(ls $DIR/fastp/*.fastq.gz)
arr=($files)
sample=${arr[$tid]}

# Run FastQC command
fastqc $sample -d . -o .

# Move all output to output folder
mv * $DIR/fastqc-trimmed


## multiqc of trimmed reads ##

#!/bin/bash
#PBS -l select=1:mem=15gb:ncpus=1
#PBS -l walltime=02:00:00
#PBS -N multiqc_trimmed_<batch>

# Set working directory
DIR=/path/to/ephemeral/CT_<batch>

# Load required module
module load anaconda3/personal
source activate cutnTag

# Run MultiQC command to match fastqc and fastqc-trimmed files 
multiqc $DIR/fastqc-trimmed/*

# Move all output to output folder
mv * $DIR/fastqc-trimmed/multiqc

```
