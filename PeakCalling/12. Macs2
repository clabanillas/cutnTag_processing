## macs2 peak calling for CUT&Tag ##

#!/bin/bash

# Set working directories
DIR=/path/to/ephemeral/CT_<batch>

# Create a separate directory inside /peakCalling for macs2 output
mkdir -p $DIR/peakCalling/macs2/unnorm_noIg/logs
mkdir -p $DIR/peakCalling/macs2/unnorm_Ig/logs



# macs2 callpeak: Unnormalized without IgG control
module load anaconda3/personal
source activate macs2

cd $DIR/peakCalling/macs2
for file in $DIR/alignment/bam/*_bowtie2.rmDup.mapped.bam
do
  sampleID=$(basename $file _bowtie2.rmDup.mapped.bam)
  OUTDIR=$DIR/peakCalling/macs2/unnorm_noIg # Define the output directory
  log=$OUTDIR/logs

  # Run macs2 callpeak with output redirection to log file in the specified output directory
  # Optional parameters from cutnrun nextflow
  macs2 callpeak -t $file -f BAMPE --keep-dup all --outdir $OUTDIR 2> "${log}/${sampleID}_macs2.log" -n "$sampleID"_nomodel --nomodel --shift -75 --extsize 150 -q 0.01
done



# macs2 callpeak: Unnormalized with IgG control for WT samples
module load anaconda3/personal
source activate macs2

cd $DIR/peakCalling/macs2
for file in $DIR/alignment/bam/*WT*_bowtie2.rmDup.mapped.bam
do
  sampleID=$(basename $file _bowtie2.rmDup.mapped.bam)
  OUTDIR=$DIR/peakCalling/macs2/unnorm_Ig # Define the output directory
  log=$OUTDIR/logs
  IgG=$DIR/alignment/bam/CT_<batch>_WT_IgG_A_bowtie2.rmDup.mapped.bam

  # Run macs2 callpeak with IgG control and output redirection to log file
  macs2 callpeak -t $file -c $IgG -f BAMPE --keep-dup all --outdir $OUTDIR 2> "${log}/${sampleID}_macs2.log" -n "$sampleID"_nextflow_Ig --nomodel --shift -75 --extsize 150 -q 0.01
done
