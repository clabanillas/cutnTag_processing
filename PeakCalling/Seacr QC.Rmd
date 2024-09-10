# Note: sample names used as examples


## 5. Assess SEACR Peak Calling quality 

```{Set up sample list paths}
## Path to the project and histone list
projPath = "/path/to/ephemeral/CT_<batch>"
conditionList = c("WT", "MSH6KO")
sampleList = c( "k27me3_A", "k27me3_B", "k36me3A_A", "k36me3A_B", "k36me3E_A", "k36me3E_B")
histList = c("k36me3", "k27me3")
```
# 5.1 Number of peaks called

```{Collect number of peaks called}
## Path to the project and histone list
projPath = "/path/to/ephemeral/CT_<batch>"
conditionList <- c("WT", "MSH6KO")
sampleList <- c("k27me3_A", "k27me3_B", "k36me3A_A", "k36me3A_B", "k36me3E_A", "k36me3E_B")
peakTypeList <- c("non_control", "top0.01", "non_relaxed")

## Initialize data frames to store peak number and width information
peakN <- data.frame()
peakWidth <- data.frame()

## Loop through conditions, histones, and peak types
for (condition in conditionList) {
  for (hist in sampleList) {
    histInfo <- strsplit(hist, "_")[[1]]
    if (histInfo[1] != "IgG") {
      for (type in peakTypeList) {
        # Construct the possible filenames
        filename_stringent <- paste0(projPath, "/peakCalling/normalised/", 
                           "CT_023_", condition, "_", hist, 
                           "_rmDup_normalised_seacr_", type, 
                           ".peaks.stringent.bed")
        
        filename_relaxed <- paste0(projPath, "/peakCalling/normalised/", 
                           "CT_023_", condition, "_", hist, 
                           "_rmDup_normalised_seacr_", type, 
                           ".peaks.relaxed.bed")
        
        # Determine which file exists
        if (file.exists(filename_stringent)) {
          filename <- filename_stringent
        } else if (file.exists(filename_relaxed)) {
          filename <- filename_relaxed
        } else {
          message(paste("No valid file exists for:", condition, hist, type))
          next
        }
        
        # Read the peak file
        peakInfo <- read.table(filename, header = FALSE, fill = TRUE) %>%
          mutate(width = abs(V3 - V2))
        
        # Store the peak number and peak width information
        peakN <- rbind(peakN, 
                       data.frame(peakN = nrow(peakInfo), 
                                  peakType = type, 
                                  Histone = histInfo[1], 
                                  Replicate = histInfo[2],
                                  Condition = condition))
        
        peakWidth <- rbind(peakWidth, 
                           data.frame(width = peakInfo$width, 
                                      peakType = type, 
                                      Histone = histInfo[1], 
                                      Replicate = histInfo[2],
                                      Condition = condition))
      }
    }
  }
}

## Save the peak number information to a CSV file
df <- data.frame(peakN %>% select(Histone, Condition, Replicate, peakType, peakN))
write.csv(df, paste("/path/to/home/", "CT_<batch>_nPeakCalls.csv", sep=""), row.names = FALSE)
```





# 5.2 Reproducibility of the peak across biological replicates

```{Calculate overlap between technical reps}
# Define the project path, condition list, histone list, and replicate list
projPath = "/rds/general/user/css119/ephemeral/CT_023"
conditionList <- c("WT", "MSH6KO")
histL <- c("k27me3", "k36me3A", "k36me3E")
repL <- c("A", "B")
peakType <- c("non_control", "top0.01", "non_relaxed")

# Initialize an empty data frame to store the results
peakOverlap <- data.frame()

# Loop through conditions, histones, peak types, and replicates
for (condition in conditionList) {
  for (type in peakType) {
    for (hist in histL) {
      overlap.gr <- GRanges()
      for (rep in repL) {
        # Construct the possible filenames
        filename_stringent <- paste0(projPath, "/peakCalling/normalised/CT_<batch>_", 
                                     condition, "_", hist, "_", rep, 
                                     "_rmDup_normalised_seacr_", type, 
                                     ".peaks.stringent.bed")
        
        filename_relaxed <- paste0(projPath, "/peakCalling/normalised/CT_<batch>_", 
                                   condition, "_", hist, "_", rep, 
                                   "_rmDup_normalised_seacr_", type, 
                                   ".peaks.relaxed.bed")
        
        # Determine which file exists
        if (file.exists(filename_stringent)) {
          filename <- filename_stringent
        } else if (file.exists(filename_relaxed)) {
          filename <- filename_relaxed
        } else {
          message(paste("No valid file exists for:", condition, hist, rep, type))
          next
        }
        
        # Read the peak file
        peakInfo <- read.table(filename, header = FALSE, fill = TRUE)
        peakInfo.gr <- GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
        
        # Find overlaps between replicates
        if(length(overlap.gr) > 0) {
          overlap.gr <- overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
        } else {
          overlap.gr <- peakInfo.gr
        }
      }
      # Use unique() to remove duplicate peak ranges
      overlap.gr = unique(overlap.gr)
      # Store the overlap information
      peakOverlap <- rbind(peakOverlap, 
                           data.frame(Condition = condition,
                                      Histone = hist, 
                                      peakType = type,
                                      peakReprod = length(overlap.gr)))
    }
  }
}

# Combine with peakN data (assuming peakN is previously generated)
peakReprod <- left_join(peakN, peakOverlap, by = c("Condition", "Histone", "peakType")) %>%
  mutate(peakReprodRate = peakReprod / peakN * 100)

# Select and save the final data
df <- peakReprod %>% 
  select(Histone, Condition, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)



write.csv(df, paste("/path/to/home/", "CT_<batch>_nPeakReprod.csv", sep=""), row.names = FALSE)
```


# 5.3 FRagment proportion in Peaks regions (FRiPs)

```{Calculate FRiP scores}

#only for _non_control_peaks#


library(chromVAR)
#paths needed
projPath = "/path/to/ephemeral/CT_<batch>"
conditionList <- c("WT", "MSH6KO")
histL <- c("k27me3", "k36me3A", "k36me3E")
repL = paste0(c("A","B"))



bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
## overlap with bam file to get count
for(condition in conditionList){
  for(hist in histL){
    for(rep in repL){
      peakRes = read.table(paste0(projPath, "/peakCalling/normalised/high_igG/", "CT_023_", condition, "_", hist, "_", rep, "_rmDup_normalised_seacr_non_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
      peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
      bamFile = paste0(bamDir, "/" , "CT_023_", condition, "_", hist, "_", rep, "_bowtie2.rmDup.mapped.bam")
      fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
      inPeakN = counts(fragment_counts)[,1] %>% sum
      inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Condition = condition, Histone = hist, Replicate = rep))
    }
  }
}


frip = left_join(inPeakData, alignResult, by = c( "Condition","Histone","Replicate")) %>% mutate(frip = inPeakN/MappedFragNum_T2T * 100)
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_T2T, AlignmentRate_T2T, FragInPeakNum = inPeakN, FRiPs = frip)

df <- data.frame(frip %>% select(Histone, Condition, Replicate, SequencingDepth, MappedFragNum_T2T, AlignmentRate_T2T, FragInPeakNum = inPeakN, FRiPs = frip))
write.csv(df, paste("/path/to/home/", "CT_<batch>_FRiPs.csv", sep=""), row.names = FALSE)
```


# 5.4 Visualisation

```{Generate plots}
fig7A = peakN %>% ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

print(fig7A)

fig7B = peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

print(fig7B)

fig7C = peakReprod %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")

print(fig7C)

fig7D = frip %>% ggplot(aes(x = Histone, y = frip, fill = Histone, label = round(frip, 2))) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks") +
  xlab("")

print(fig7D)

ggarrange(fig7A, fig7B, fig7C, fig7D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
ggsave(paste0("/path/to/home/", "CT_<batch>_peakCalling_stats.pdf"), height = 10, width = 20)
```




