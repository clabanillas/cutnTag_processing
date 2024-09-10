---
title: "Scale Factor of <batch> Samples"
author: "Your Name"
date: "`r Sys.Date()`"
output: html_document
---

```{Load libraries}
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR)  # For FRiP analysis and differential analysis
library(DESeq2)    # For differential analysis section
library(ggpubr)    # For customizing figures
library(corrplot)  # For correlation plot
```

## Sample and Histone List Setup

```{Sample list set up}
projPath <- "/path/to/ephemeral/CT_<batch>"
sampleList <- c(
  "CT_023_MSH6KO_IgG_A", "CT_023_MSH6KO_IgG_B", "CT_023_MSH6KO_k27me3_A",  "CT_023_MSH6KO_k27me3_B", 
  "CT_023_MSH6KO_k36me3A_A", "CT_023_MSH6KO_k36me3A_B", "CT_023_MSH6KO_k36me3E_A", "CT_023_MSH6KO_k36me3E_B", 
  "CT_023_WT_IgG_A", "CT_023_WT_IgG_B", "CT_023_WT_k27me3_A", "CT_023_WT_k27me3_B", 
  "CT_023_WT_k36me3A_A", "CT_023_WT_k36me3A_B", "CT_023_WT_k36me3E_A", "CT_023_WT_k36me3E_B"
)
histList <- c("IgG", "k36me3A", "k36me3E", "k27me3")
```

## Collecting Alignment Results from Bowtie2 Summary Files

```{Collect alignment results in dataframe}
alignResult <- data.frame()
for(hist in sampleList) {
  alignRes <- read.table(
    paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2.txt"), 
    header = FALSE, fill = TRUE
  )
  
  alignRate <- substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo <- strsplit(hist, "_")[[1]]
  
  alignResult <- data.frame(
    Sample = paste(histInfo[1], histInfo[2], histInfo[3], histInfo[4], histInfo[5], sep = "_"), 
    Condition = histInfo[3], 
    Histone = histInfo[4], 
    Replicate = histInfo[5], 
    SequencingDepth = as.numeric(alignRes$V1[1]), 
    MappedFragNum_T2T = as.numeric(alignRes$V1[4]) + as.numeric(alignRes$V1[5]), 
    AlignmentRate_T2T = as.numeric(alignRate)
  ) %>% rbind(alignResult, .)
}
alignResult$Histone <- factor(alignResult$Histone, levels = histList)
alignResult <- alignResult %>% mutate(AlignmentRate_T2T = paste0(AlignmentRate_T2T, "%"))
alignResult
```

## Summarizing Duplication Information from Picard Summary Outputs

```{Collect duplication results an djoin with alignment results}
dupResult <- data.frame()
for(hist in sampleList) {
  dupRes <- read.table(
    paste0(projPath, "/alignment/removeDuplicate/picard_summary/", hist, "_picard.rmDup.txt"), 
    header = TRUE, fill = TRUE
  )
  histInfo <- strsplit(hist, "_")[[1]]
  
  dupResult <- data.frame(
    Sample = paste(histInfo[1], histInfo[2], histInfo[3], histInfo[4], histInfo[5], sep = "_"), 
    Condition = histInfo[3], 
    Histone = histInfo[4], 
    Replicate = histInfo[5], 
    MappedFragNum_T2T = as.numeric(dupRes$READ_PAIRS_EXAMINED[1]), 
    DuplicationRate = as.numeric(dupRes$PERCENT_DUPLICATION[1]) * 100, 
    EstimatedLibrarySize = as.numeric(dupRes$ESTIMATED_LIBRARY_SIZE[1])
  ) %>% mutate(UniqueFragNum = MappedFragNum_T2T * (1 - DuplicationRate / 100)) %>% rbind(dupResult, .)
}
dupResult$Histone <- factor(dupResult$Histone, levels = histList)
alignDupSummary <- left_join(alignResult, dupResult, by = c("Sample", "Condition", "Histone", "Replicate", "MappedFragNum_T2T")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))
alignDupSummary
```

## Processing Spike-in counts 

```{Process barcode counts and calculate scale factor}
# Read barcode counts
barcodes_txt <- read.delim("/path/to/ephemeral/CT_<batch>/kmet/kmet_barcode_counts_CT_<batch>.txt")
spike_in_reads <- read.csv("/path/to/ephemeral/CT_<batch>/kmet/kmet_barcode_counts_CT_<batch>.txt", sep = "\t")

# Sum barcode counts for paired reads
spike_in_reads_pair <- spike_in_reads[, c(T, F)] + spike_in_reads[, c(F, T)]

# Process row names
split_names <- strsplit(rownames(spike_in_reads), split = "_")
ID_files <- sapply(split_names, function(x) paste(x[-length(x)], collapse = "_"))
ID_files <- ID_files[c(T, F)]
rownames(spike_in_reads_pair) <- ID_files

# IMPORTANT: Match with alignDupSummary rows
# Set row names for alignDupsummary based on a column
row.names(alignDupSummary) <- alignDupSummary$Sample
#row.names(spike_in_reads_pair) already set and are the same
# Order one dataframe by the row names of the other (alignDupSummary)
spike_in_reads_pair <- spike_in_reads_pair[match(row.names(alignDupSummary), row.names(spike_in_reads_pair)), ]


# Generate spike-in dataframe
spike_in_reads_pair <- as.data.frame(spike_in_reads_pair)
spike_in_reads_pair[,"sum"] <- rowSums(spike_in_reads_pair)
spike_in_reads_pair[,"Histone"] <- alignDupSummary$Histone
spike_in_reads_pair[,"Technical_Replicate"] <- alignDupSummary$Replicate
spike_in_reads_pair[,"Barcode_Interest"] <- apply(spike_in_reads_pair, 1, FUN = function(x) {
  hist <- x["Histone"]
  result <- x[grepl(hist, names(x))]
  return(as.numeric(result))
})
spike_in_reads_pair[,"Total_Reads"] <- alignDupSummary$SequencingDepth
spike_in_reads_pair[,"Percentage_Spike_In"] <- spike_in_reads_pair$sum / spike_in_reads_pair$Total_Reads * 100

# Calculate Scale Factor and Specificity
spike_in_reads_pair[,"Spike_In_Factor"] <- 1 / spike_in_reads_pair$sum
spike_in_reads_pair[,"Specificity_Factor"] <- spike_in_reads_pair$Barcode_Interest / spike_in_reads_pair$sum
spike_in_reads_pair[,"Final_Factor"] <- spike_in_reads_pair$Spike_In_Factor * spike_in_reads_pair$Specificity_Factor * 1e6

# Output scale factors
write.csv(spike_in_reads_pair$Final_Factor, "/path/to/ephemeral/CT_<batch>/kmet/scale_factors_CT_<batch>.csv", row.names = FALSE)

# Output processed information
write.csv(spike_in_reads_pair, "/path/to/ephemeral/CT_<batch>/kmet/Proc_Info_CT_<batch>.csv", row.names = FALSE)
```




