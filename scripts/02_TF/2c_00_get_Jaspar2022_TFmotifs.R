rm(list = ls())
setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/utils/") # original: ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/

library(readxl)

#BiocManager::install("MotifDb")
library(MotifDb)
MotifDb_names = names(MotifDb)
idx_JASPAR2022 = which(grepl("Hsapiens-jaspar2022", MotifDb_names))
JASPAR2022_names = MotifDb_names[idx_JASPAR2022] # 691
mf = MotifDb[JASPAR2022_names]
export(mf, "Jaspar2022motifs.jaspar", 'jaspar')
export(mf, "Jaspar2022motifs.meme", 'meme')

