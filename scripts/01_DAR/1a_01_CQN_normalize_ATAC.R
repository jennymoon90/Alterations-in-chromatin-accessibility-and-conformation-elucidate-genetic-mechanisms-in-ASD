### Goal: Use CQN to calculate normalization factors for each peak accounting for GC content, peak width, and total number of unique non-mitochondrial fragments sequenced. 

rm(list = ls())
options(stringsAsFactors = F)
# BiocManager::install("cqn")
library(cqn)
library(tidyverse)
# library(scales)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")

########## pipeline ########## 
### 1. Load data
  ## 1) ATAC peak raw counts
  ## 2) ATAC peak GC contents
  ## 3) ATAC library biological and technical covariates
### 2. Run CQN
### 3. MA-plot (Figure 1c)
##############################

### 1. Load data
## 1) ATAC peak raw counts:
lnames = load("../../data/ATAC_raw_counts.rda") # df, original file name 'Peaksets_minOverlap5_chr1to22XY_featureCounts.rda'
counts_tbl = df; rm(df)
rownames(counts_tbl) = paste0(counts_tbl$Chr, "_", counts_tbl$Start, "_", counts_tbl$End)

## 2) ATAC peak GC contents:
# Run 1a_02_bedtools_nuc_GCcontent.sh to get GC content of the ATAC peaks. Load the resulting .tsv file:
peak_gc = read.table("ATAC_bedtoolsnuc.tsv", header = T) # original file name 'Peaksets_minOverlap5_bedtoolsnuc.tsv'
colnames(peak_gc)[1:5] = c("Chr", "Start", "End", "pct_AT", "pct_GC")
peak_info = left_join(counts_tbl[,1:6], peak_gc[,1:5])
rownames(peak_info) = rownames(counts_tbl)

peak_info = peak_info[,c(6,8)] # 128036 peaks x 2 columns (Length and pct_GC)
counts_tbl = counts_tbl[,7:ncol(counts_tbl)] # 128036 peaks x 38 samples

## 3) ATAC library biological and technical covariates:
lnames = load("../../data/ATAC_library_covariates.rda") # df_wAllCov. original file name: 3_QC/BiolTechCov_and_QC.rda
Covariates = df_wAllCov; rm(df_wAllCov)
Covariates = Covariates[match(colnames(counts_tbl), Covariates$sample_name),] # 38 rows/samples

nonMitoUniqRd = Covariates$N_nonmito_unique_reads
names(nonMitoUniqRd) = Covariates$sample_name

### 2. Run CQN
stopifnot(all(rownames(counts_tbl) == rownames(peak_info)))
stopifnot(colnames(counts_tbl) == names(nonMitoUniqRd))

cqn.fit <- cqn(counts_tbl, lengths = peak_info$Length,
                  x = peak_info$pct_GC, sizeFactors = nonMitoUniqRd,
                  verbose = TRUE)
# 2022/08/12 00:28-00:29 AM

## Examine CQN plots
pdf("../../plots/01_DAR/1a_01_CQN.pdf", width = 8, height = 5)
par(mfrow=c(1,2))
cqnplot(cqn.fit, n = 1, xlab = "GC content", lty = 1, ylim = c(-1.5,5)) # 
cqnplot(cqn.fit, n = 2, xlab = "length", lty = 1, ylim = c(-1.5,5)) # 
dev.off()

## Normalized peak counts
RPKM.cqn <- cqn.fit$y + cqn.fit$offset # These values are on the log2-scale
save(RPKM.cqn, cqn.fit, Covariates, peak_info, counts_tbl,
     file = "1a_01_CQN.rda")


# ------------ the end ------

### 3. MA-plot (Figure 1c)
# ASD differential ATAC peaks with FDR < 0.05
rm(list = ls())
TextSize = 5 # 5 pt, eg. text = element_text(size = TextSize)
AnnoSize = 1.8 # font 0.093 inch tall or 5 pt, used in annotate(size = AnnoSize) or geom_text(size = AnnoSize)
PointSize = 1 # used as PointSize/2 or PointSize/3 in geom_point

lnames = load("1a_01_CQN.rda") # "RPKM.cqn"   "cqn.fit"    "Covariates" "peak_info"  "counts_tbl"
lnames = load(paste0(DIR,"5_DiffATAC/DiffATAC.rda")) # DiffATAC
stopifnot(rownames(RPKM.cqn) == rownames(DiffATAC))
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

plot(MA_df$Mean, MA_df$ATAC_logFC, pch = 19, col = MA_df$col)
# idenify out skirts for idx_keep
abline(h = 0.65)
abline(h = -0.65)
abline(v = -2)

# to reduce pdf size
idx_keep = which(MA_df$ATAC_logFC >= 0.65 | MA_df$ATAC_logFC <= -0.65 | MA_df$Mean <= -2 | MA_df$Mean >= 3) # 177 data points
idx1 = which(MA_df$Mean >= 2) # 2424
idx2 = which(MA_df$Mean < 2.5) # 127396
idx2 = idx2[which(idx2 <= tmp$Freq[tmp$Var1 == "lightgrey"])] # 122392 n.s. data points with x-axis < 2.5
set.seed(42)
idx2s = sample(idx2, 10e3) # tried seq(5e3, 30e3, by = 5e3), 10e3 is sufficient for 2.3 x 2 inches of plot, reduced pdf file size from 6.8MB to 1.1MB
MA_df_sampled = MA_df[c(idx_keep, idx1, idx2s, (tmp$Freq[tmp$Var1 == "lightgrey"]+1):nrow(MA_df)),]

Fig1c = 
  # MA_df %>%
  MA_df_sampled %>% # reduced plot size
  ggplot(aes(x = Mean, y = ATAC_logFC, col = col)) +
  geom_point(size = PointSize/2) +
  scale_color_manual(values = c("red" = "red", "blue" = "blue", "lightgrey" = "lightgrey")) + 
  annotate("text", x = 2, y = 1.1, label = paste0("Up-reg: ", tmp$Freq[tmp$Var1 == "red"], " peaks\nDown-reg: ", tmp$Freq[tmp$Var1 == "blue"], " peaks\n(FDR < ",fdr_t,")"), size = AnnoSize) +
  theme_bw() + # size = 2 gives 8 pt
  #labs(col = "DARs") +
  xlab("Log2 mean accessibility") +
  ylab("Log2 fold change") +
  ggtitle(paste0("Total ", nrow(MA_df), " ATAC peaks")) +
  theme(legend.position = "none", #c(0.9,0.8),
        text = element_text(size = TextSize),
        axis.text = element_text(size = TextSize), # required for axis text to be 5 pt
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.5, size = TextSize),
        # element_blank() removes grey grids but also axis ticks
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

pdf("1a_01_DAR_MA_plot_Fig1c.pdf", width = 2.3, height = 2) # original file name: Figure1c_ATAC_MAplot.pdf
print(Fig1c)
dev.off()

pdf("")