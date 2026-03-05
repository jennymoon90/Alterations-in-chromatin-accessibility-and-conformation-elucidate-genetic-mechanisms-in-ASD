rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(Repitools)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")

########## pipeline ########## 
### 1. Load data
  ## 1) Load differential ATAC results
  ## 2) Load differential H3K27Ac results
### 2. Correlation between logFC of differential ATAC and H3K27ac peaks
##############################

### 1. Load data
## 1) Load differential ATAC results
load("1a_03_DiffATAC.rda") # DiffATAC
DiffATAC_location = DiffATAC[,1:3]

## 2) Load differential H3K27Ac results
lnames = load("../../data/DiffH3K27ac_Ramaswami2020NatCommun.Rdata") # DAR. Original file: /Users/dhglab/Documents/Documents/Geschwind_lab/LAB/Database_download/32_Gokul_H3K27Ac_datasets/processed_data/differential_acetylation/Convergent.Rdata

DAR_location = as.data.frame(str_split_fixed(rownames(DAR), "-", 3))
colnames(DAR_location) = c("chr", "start", "end")
DAR_location$start = as.integer(DAR_location$start)
DAR_location$end = as.integer(DAR_location$end)
DAR = cbind(DAR_location, DAR)
rownames(DAR_location) = rownames(DAR)

unique(DiffATAC$chr) # chr1-Y
unique(DAR$chr) # chr1-Y, chrM and xxx_random
DAR = DAR[-which(grepl("_random", DAR$chr)),] # 56451 H3K27ac peaks
DAR = DAR[DAR$chr != "chrM",] # 56449 H3K27ac peaks

### 2. Correlation between logFC of differential ATAC and H3K27ac peaks
# Filter for FDR < 0.05 of all profiles
DiffATAC = DiffATAC[!is.na(DiffATAC$ATAC_FDR),] # 127971 ATAC peaks
DAR = DAR[!is.na(DAR$beta.condition),] # 56168 H3K27ac peaks
sigATAC = DiffATAC[DiffATAC$ATAC_FDR < 0.05,] # 5033 significant ATAC peaks
sigH3K27ac = DAR[DAR$p.condition.fdr < 0.05,] # 8158 significant H3K27ac peaks

# Intersect the two significant peak sets
sigATAC_gr = annoDF2GR(sigATAC)
sigH3K27ac_gr = annoDF2GR(sigH3K27ac)
hits = as.data.frame(findOverlaps(sigATAC_gr, sigH3K27ac_gr))
hits$ATACpeak = rownames(sigATAC)[hits$queryHits]
hits$H3K27ACpeak = rownames(sigH3K27ac)[hits$subjectHits]
hits$ATAC_logFC = sigATAC$ATAC_logFC[hits$queryHits]
hits$H3K27ac_logFC = sigH3K27ac$beta.condition[hits$subjectHits]

# Calculate the correlation between ATAC and H3K27ac logFCs
hit_cor = round(cor(hits$ATAC_logFC, hits$H3K27ac_logFC),2) # 0.82
tmp = cor.test(hits$ATAC_logFC, hits$H3K27ac_logFC)
hit_corsig = formatC(tmp$p.value, 1) # 4e-197

df = data_frame(x = hits$ATAC_logFC, y = hits$H3K27ac_logFC)
model = lm(y ~ 0 + x, df)
model_res = summary(model)
r_coef = round(sqrt(model_res$r.squared), digits = 2) # 0.89
p_coef = formatC(model_res$coefficients[1,4], digits = 1) # 1e-275

# Plot the correlation
TextSize = 5
AnnoSize = 1.8
PointSize = 1

Fig1d = 
  hits %>%
  ggplot(aes(x = ATAC_logFC, y = H3K27ac_logFC)) +
  geom_point(size = PointSize/2) +
  geom_abline(intercept = 0, slope = lm_coef, col = "red") +
  annotate("text", x = 0.4, y = -0.05, label = paste0("R^2 = ", round(r_coef^2, 1), "\np = ", p_coef), col = "red", size = AnnoSize) +
  theme_bw() +
  xlab("ATAC peak logFC") +
  ylab("H3K27ac peak logFC") +
  ggtitle(paste0(length(unique(hits$ATACpeak)), " diff. ATAC overlap with ", length(unique(hits$H3K27ACpeak)), " diff. H3K27ac peaks")) +
  theme(text = element_text(size = TextSize),
        axis.text = element_text(size = TextSize),
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.5, size = TextSize),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

pdf("../../plots/01_DAR/1b_01_cor_logFC_DiffATAC_H3K27ac_Fig1d.pdf", width = 2.3, height = 2) # original file name: Figure1d_cor_diffATAC_diffH3K27Ac.pdf
print(Fig1d)
dev.off()

