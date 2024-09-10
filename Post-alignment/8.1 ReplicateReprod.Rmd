# Note: Refer to script 5.1 for previous steps 0-2 and 6.1 for previous step 3


## 4. Assessment of replicate reproducibility 

```{Collect bin500 files and build count matrix}
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
reprod = c()
fragCount = NULL
for(hist in sampleList){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(projPath, "/alignment/bed/", hist, "_rmDup_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
    
  }else{
    
    fragCountTmp = read.table(paste0(projPath, "/alignment/bed/", hist, "_rmDup_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "pairwise.complete.obs", method = "spearman") 

pdf(paste0("/path/to/home/", filename, "_repcorrplot.pdf"), width = 12, height = 10)

fig6 <- corrplot(M, 
           method = "color", 
           outline = T, 
           addgrid.col = "darkgray", 
           order="hclust", addrect = 3, 
           rect.col = "black", 
           rect.lwd = 3,
           cl.pos = "b", 
           tl.col = "indianred4", 
           tl.cex = 1, 
           cl.cex = 1, 
           addCoef.col = "black", 
           number.digits = 2, 
           number.cex = 1, 
           col = colorRampPalette(c("midnightblue","white","darkred"))(100))

dev.off()
```

