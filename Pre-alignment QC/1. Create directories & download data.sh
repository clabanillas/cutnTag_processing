```bash

# Note: Fastq file names were changed to the EPIMEM convention CT_<batch>_<day>_<histone>_<rep>

# Define base directories with placeholders
ephem=/path/to/ephemeral
epimem=/path/to/project/epimem
data=/path/to/project/epimem/cutntag

# Define specific directories with placeholders
DIR=$ephem/DIR
scripts=$DIR/scripts

# Make directories
mkdir -p $DIR/raw
mkdir -p $DIR/fastqc
mkdir -p $DIR/kmet
mkdir -p $DIR/fastqc/multiqc
mkdir -p $DIR/fastp
mkdir -p $DIR/fastp/multiqc
mkdir -p $DIR/fastqc-trimmed
mkdir -p $DIR/fastqc-trimmed/multiqc

# Alignment directories
mkdir -p $DIR/alignment/sam/bowtie2_summary
mkdir -p $DIR/alignment/bam
mkdir -p $DIR/alignment/bed
mkdir -p $DIR/alignment/bedgraph
mkdir -p $DIR/alignment/bedgraph/histControl
mkdir -p $DIR/alignment/bedgraph/normalised
mkdir -p $DIR/alignment/bedgraph/normalised/histControl
mkdir -p $DIR/alignment/removeDuplicate/picard_summary
mkdir -p $DIR/alignment/sam/fragmentLen
mkdir -p $DIR/alignment/bigwig

# Peak directories
mkdir -p $DIR/peakCalling
mkdir -p $DIR/peakCalling/normalised
mkdir -p $DIR/peakCalling/bigwig
mkdir -p $DIR/peakCalling/normalised/bigwig

## Getting the data ##

# Copy over data 
cp -r $data/raw $DIR

# Check data integrity
md5sum -c MD5.txt 

# Uncompress data
tar -xvf <file>.tar

## Changing names of files to EPIMEM nomenclature for clarity ##
# Format: <batchID>_<biorep>_<mark>_<rep>
# Example: <CT_010>_<D43B>_<k27me3>_<A>_
```
