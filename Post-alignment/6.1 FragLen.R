# Note: Refer to script 5.1 for previous steps 0-2

## 3. Assessing fragment length size distribution 

```{Collect the fragment size information}
fragLen = c()
for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[4], Replicate = histInfo[5], sampleInfo = hist) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)
```

```{Generate plots}
##PLOTS
## Generate the fragment size density plot (violin plot)

fig5A = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  theme(plot.margin = margin(t = 1,  # Top margin
                             r = 1,  # Right margin
                             b = 3,  # Bottom margin
                             l = 5,  # Left margin
                             unit = "cm")) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

print(fig5A)

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

print(fig5B)
```

```{Combine plots to PDF}
ggarrange(fig5A, fig5B, ncol = 2)
ggsave(paste0("/path/to/home/", filename, "_fragLen_stats.pdf"), height = 10, width = 20)
```
