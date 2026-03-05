### Goal: MA plot of the differential ATAC results by Diagnosis (ASD vs. CTL) 

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)

TextSize = 5 # 5 pt, eg. text = element_text(size = TextSize)
AnnoSize = 1.8 # font 0.093 inch tall or 5 pt, used in annotate(size = AnnoSize) or geom_text(size = AnnoSize)
PointSize = 1 # used as PointSize/2 or PointSize/3 in geom_point

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")

########## pipeline ########## 
### 1. Load data
  ## 1) CQN normalized ATAC peak signal
  ## 2) ATAC differential analysis results
### 2. MA-plot (Figure 1c)
  ## 1) Full plot  
  ## 2) To reduce pdf size, randomly sample the data points
##############################

### 1. Load ATAC differential analysis results
## 1) CQN normalized ATAC peak signal
lnames = load("1a_01_CQN.rda") # "RPKM.cqn"   "cqn.fit"    "Covariates" "peak_info"  "counts_tbl"
## 2) ATAC differential analysis results
lnames = load("1a_03_DiffATAC.rda") # DiffATAC
stopifnot(rownames(RPKM.cqn) == rownames(DiffATAC))

### 2. MA-plot (Figure 1c)
# Highlight ASD differential ATAC peaks with FDR < 0.05
Mean = apply(RPKM.cqn, 1, mean)
MA_df = cbind(DiffATAC, Mean)

fdr_t = 0.05
logfc_t = 0

MA_df$col = case_when(
  MA_df$ATAC_FDR >= fdr_t | is.na(MA_df$ATAC_FDR) ~ "lightgrey", # "n.s.",
  MA_df$ATAC_FDR < fdr_t & MA_df$ATAC_logFC > logfc_t ~ "red", # "Up-reg",
  MA_df$ATAC_FDR < fdr_t & MA_df$ATAC_logFC < -logfc_t ~ "blue" # "Down-reg" 
)

tmp = as.data.frame(table(MA_df$col)) # blue (down) 2248, red (up) 2785, lightgrey 122938
MA_df = MA_df %>%
  arrange(desc(ATAC_FDR))

## 1) Full plot  
plot(MA_df$Mean, MA_df$ATAC_logFC, pch = 19, col = MA_df$col)
# idenify out skirts for idx_keep
abline(h = 0.65)
abline(h = -0.65)
abline(v = -2)

## 2) To reduce pdf size, randomly sample the data points
idx_keep = which(MA_df$ATAC_logFC >= 0.65 | MA_df$ATAC_logFC <= -0.65 | MA_df$Mean <= -2 | MA_df$Mean >= 3) # 175 data points
idx1 = which(MA_df$Mean >= 2) # 2429
idx2 = which(MA_df$Mean < 2.5) # 122396
idx2 = idx2[which(idx2 <= tmp$Freq[tmp$Var1 == "lightgrey"])] # 122392 n.s. data points with x-axis < 2.5
set.seed(42)
idx2s = sample(idx2, 10e3) # tried seq(5e3, 30e3, by = 5e3), 10e3 is sufficient for 2.3 x 2 inches of plot, reduced pdf file size from 6.8MB to 1.1MB
MA_df_sampled = MA_df[c(idx_keep, idx1, idx2s, (tmp$Freq[tmp$Var1 == "lightgrey"]+1):nrow(MA_df)),]

Fig1c = 
  MA_df_sampled %>% # reduced plot size
  ggplot(aes(x = Mean, y = ATAC_logFC, col = col)) +
  geom_point(size = PointSize/2) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "lightgrey" = "lightgrey")) + 
  annotate("text", x = 2, y = 1.1, label = paste0("Up-reg: ", tmp$Freq[tmp$Var1 == "red"], " peaks\nDown-reg: ", tmp$Freq[tmp$Var1 == "blue"], " peaks\n(FDR < ",fdr_t,")"), size = AnnoSize) +
  theme_bw() + 
  xlab("Log2 mean accessibility") +
  ylab("Log2 fold change") +
  ggtitle(paste0("Total ", nrow(MA_df), " ATAC peaks")) +
  theme(legend.position = "none", 
        text = element_text(size = TextSize),
        axis.text = element_text(size = TextSize), # required for axis text to be 5 pt
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.5, size = TextSize),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

pdf("../../plots/01_DAR/1a_04_DAR_MAplot_Fig1c.pdf", width = 2.3, height = 2) # original file name: Figure1c_ATAC_MAplot.pdf
print(Fig1c)
dev.off()
