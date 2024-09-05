```bash

## kmet loop for barcode counts ##

# kmetloop.sh is a script that counts barcodes corresponding to the Epycypher k-methyl modification panel spike-in from fastq files within a folder
# Remember to customize PBS -N and PBS -J 0-N 

#!/bin/bash
#PBS -l select=1:mem=15gb:ncpus=1
#PBS -l walltime=02:00:00
#PBS -N kmetloop_<batch>
#PBS -J 0-15

tid=$PBS_ARRAY_INDEX

# Set working directory containing all the fq.gz files
DIR=/path/to/ephemeral/CT_<batch>
OUTDIR=$DIR/kmet

# List your input files
files=$(ls $DIR/*/*/*/*/*_1.fq.gz)
arr=($files)
sample=${arr[$tid]}

# Create sampleID, common for matched _1 and _2 files 
sampleID=$(basename $sample _1.fq.gz)

# Unzip fastq files first
gunzip -k $DIR/*/*/*/*/"$sampleID"_1.fq.gz
gunzip -k $DIR/*/*/*/*/"$sampleID"_2.fq.gz

# Run kmetloop.sh which counts barcodes sequences corresponding to kmet panel spike-in from sequencing data

read1=$DIR/*/*/*/*/"$sampleID"_1.fq
read2=$DIR/*/*/*/*/"$sampleID"_2.fq
samplename=$(basename $read1 _1.fq)

echo $(basename $read1) > "$samplename"_R1.txt
for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA; do
    grep -c $barcode $read1 >> "$samplename"_R1.txt  
done

echo $(basename $read2) >> "$samplename"_R2.txt
for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA; do
    grep -c $barcode $read2 >> "$samplename"_R2.txt
done 

# Paste both reads into a single file
paste "$samplename"_R1.txt "$samplename"_R2.txt | column -s '' -t > "$samplename"_bc_counts.txt

# Move all output files to output folder
mv *_bc_counts.txt $OUTDIR


## TERMINAL ##
# Collate all barcode counts into a single .txt file
cd /path/to/ephemeral/CT_<batch>/kmet
LABEL=kmet_barcode_counts_CT_<batch>_
paste *_bc_counts.txt | column -s '' -t > "$LABEL".txt

```
