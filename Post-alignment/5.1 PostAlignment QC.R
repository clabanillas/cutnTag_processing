---
title: "Post_Alignment_QC_CT_<batch>"
author: "your name"
date: "`r Sys.Date()`"
output: html_document
---

## 0. Load libraries and R working environment 

```{Load libraries}
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(corrplot) ## For correlation plot
```


```{Set up sample list}
## Path to the project and histone list
projPath = "/path/to/ephemeral/CT_<batch>"
sampleList=c("CT_023_MSH6KO_IgG_A", "CT_023_MSH6KO_IgG_B", "CT_023_MSH6KO_k27me3_A",  "CT_023_MSH6KO_k27me3_B", "CT_023_MSH6KO_k36me3A_A", "CT_023_MSH6KO_k36me3A_B", "CT_023_MSH6KO_k36me3E_A", "CT_023_MSH6KO_k36me3E_B", "CT_023_WT_IgG_A", "CT_023_WT_IgG_B", "CT_023_WT_k27me3_A", "CT_023_WT_k27me3_B", "CT_023_WT_k36me3A_A", "CT_023_WT_k36me3A_B", "CT_023_WT_k36me3E_A", "CT_023_WT_k36me3E_B")
histList = c("IgG", "k36me3A", "k36me3E", "k27me3")
#assign a filename to the batch name
filename = "CT_023"
```

## 1. Assessment of Alignment quality 

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

```{Generate alignment QC plots}
#Sequencing depth plot
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Condition), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Sequencing Depth per Million") +
  xlab("") + 
  ggtitle("A. Sequencing Depth")

print(fig3A)

#Mapped fragments per million
fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_T2T/1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Condition), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Mapped Fragments per Million") +
  xlab("") +
  ggtitle("B. Alignable Fragment (T2T)")

print(fig3B)

#Alignment rate
fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_T2T, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Condition), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Mapped Fragments") +
  xlab("") +
  ggtitle("C. Alignment Rate (T2T)")

print(fig3C)
```

```{Combine plots into a single pdf}
ggarrange(fig3A, fig3B, fig3C, ncol = 3, nrow=1, common.legend = TRUE, legend="bottom")
ggsave(paste0("/path/to/home/", filename, "_alignment_stats.pdf"), height = 10, width = 20)
```

## 2. Assessment of Duplication Rate 

## Summarizing Duplication Information from Picard Summary Outputs

```{Collect duplication results and join with alignment results}
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

```{Generate duplication QC plots}
fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Condition), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Duplication Rate (*100%)") +
  xlab("") 

print(fig4A)

#estimated library size
fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Condition), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Estimated Library Size") +
  xlab("") 

print(fig4B)

#number of unique fragments
fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Condition), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("# of Unique Fragments") +
  xlab("")

print(fig4C)


```{Combine plots into a single pdf}
ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
ggsave(paste0("/path/to/home/", filename, "_duplication_stats.pdf"), height = 10, width = 20)
```
