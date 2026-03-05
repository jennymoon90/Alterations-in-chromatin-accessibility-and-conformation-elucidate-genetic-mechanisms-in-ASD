# qrsh -l h_data=10G,h_rt=8:00:00
# module load R # Loading R/4.0.2
# R

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/2c_01_TOBIAS/ind_BINDetect_Jaspar2022/all_bindetect_results/") # original: u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/ind_BINDetect_Jaspar2022/all_bindetect_results/

all_files = list.files(pattern = ".*.txt"); length(all_files) # 38

df = read.table(all_files[1], header = T)
str(df) # 691 motifs x 7 columns
colnames(df)[c(3,6)]  # "motif_id" "X2GGP.1_mean_score"
df = df[,c(3,6)]
  
for (f in all_files[2:length(all_files)]) {
  ind_df = read.table(f, header = T)
  df = cbind(df, ind_df[,6])
  colnames(df)[ncol(df)] = colnames(ind_df)[6]
}

df = as.data.frame(df)
str(df)
save(df, file = "../../../2c_02_Jaspar2022_TFmeanScore_by_sample.rda") # original: Jaspar2022hs_TFmotifs_mean_score_by_sample.rda
