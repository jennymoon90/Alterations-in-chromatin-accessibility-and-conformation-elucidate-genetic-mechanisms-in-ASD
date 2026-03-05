### Goal of this set of scripts (2e_01): Show that ATAC peaks bound by TFs with differential GE going in the same direction as binding score are more likely to be differentially accessible in the same direction. 
# 1) Extract TF-bound ATAC peaks in ASD vs CTL -> quantify N_peaks bound 
# 2) Randomly draw n ATAC peaks -> build null distribution of probability of ATAC diff in the same dir -> p-val of observed probablity of ATAC bound being diff in the same dir 
# 3) heatmap of ATAC log(CPM+5) in ASD and CTL samples (row x col = peak x sample, left: all CTL samples, right: all ASD samples)

# This script focuses on 1) Number of ATAC peaks bound by TF of interest in ASD vs CTL

########## pipeline ########## 
### 1. Load TF of interest (TFoi) and covariates
# TFoi: TFs with their expression significantly changing in the same orientation as their binding score
### 2. Get the number of ATAC peaks bound by each TFoi in each sample
### 3. Label the TF changes (up/down)
### 4. Boxplot separately for TFs with increased vs decreased binding
##############################


########################
## On Hoffman2 server ##
########################

# qrsh -l h_data=10G,h_rt=8:00:00
# module load R # Loading R/4.0.2
# R

rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(reshape2)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/") # original directory: /u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/ind_BINDetect_Jaspar2022/all_bindetect_results

### 1. Load TF of interest (TFoi) and covariates
lnames = load("../../data/ATAC_library_covariates.rda") # df_wAllCov. # original: ../../3_QC/BiolTechCov_and_QC.rda
Covariates = df_wAllCov; rm(df_wAllCov)
Covariates = Covariates[Covariates$QC_summary == "good",]
ASD_samples = Covariates$sample_name[Covariates$Diagnosis == "ASD"]
CTL_samples = Covariates$sample_name[Covariates$Diagnosis == "CTL"]

lnames = load("2d_TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda") # original file: TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda
TFoi = TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF

### 2. Get the number of ATAC peaks bound by each TFoi in each sample
df = matrix(nrow = length(TFoi), ncol = 0)
DIR_files = "2c_01_TOBIAS/ind_BINDetect_Jaspar2022/all_bindetect_results/"
for (s in c(ASD_samples, CTL_samples)) {
  f = paste0(DIR_files, s, "_bindetect_results.txt")
  res = read.table(f, header = T)
  res_TFoi = res[grepl(paste0(TFoi, collapse = "|"), res$motif_id),]
  df = cbind(df, res_TFoi[,ncol(res_TFoi)])
  rownames(df) = res_TFoi$motif_id
  colnames(df)[ncol(df)] = s
  print(df)
}
Npeaks_BoundByTF_EaSample = df; rm(df)

tmp = str_split_fixed(rownames(Npeaks_BoundByTF_EaSample), "-", 3)
rownames(Npeaks_BoundByTF_EaSample) = tmp[,3]

df_melt = melt(Npeaks_BoundByTF_EaSample)
str(df_melt)
colnames(df_melt) = c("TFmotif", "sample_name", "Npeak_bound")
df_melt$Diagnosis = Covariates$Diagnosis[match(df_melt$sample_name, Covariates$sample_name)]
Npeaks_BoundByTF_EaSample_melt = df_melt; rm(df_melt)

save(Npeaks_BoundByTF_EaSample, Npeaks_BoundByTF_EaSample_melt, file = "2e_01_Npeaks_boundByTF_EaSample.rda") # original file name: Npeaks_BoundByTF_EaSample.rda

#####################
## Now work on MAC ##
#####################

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggforce)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/") # original directory: /u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/ind_BINDetect_Jaspar2022/

### 3. Label the TF changes (up/down)
lnames = load("2d_diffTOBIAS_DGE.rda") # "DiffBS_lm_Volcano" "TFclusters". # original file: DiffBS_lm_Volcano_wCluster_DGE.rda
lnames = load("2e_01_Npeaks_boundByTF_EaSample.rda") # "Npeaks_BoundByTF_EaSample"      "Npeaks_BoundByTF_EaSample_melt". 
Npeaks_BoundByTF_EaSample_melt$change = TFclusters$change[match(Npeaks_BoundByTF_EaSample_melt$TFmotif, TFclusters$TFmotif)]
Npeaks_BoundByTF_EaSample_melt = Npeaks_BoundByTF_EaSample_melt %>%
  arrange(change, TFmotif)
Npeaks_BoundByTF_EaSample_melt$TFmotif = factor(Npeaks_BoundByTF_EaSample_melt$TFmotif, levels = unique(Npeaks_BoundByTF_EaSample_melt$TFmotif))
save(Npeaks_BoundByTF_EaSample_melt, file = "2e_01_num_ATACpeaks_boundBy_TFoi.rda") # original file name: Number_of_ATACpeaks_BoundBy_TF_BindingScoreGoingInSameDirWithDGE.rda

### 4. Boxplot separately for TFs with increased vs decreased binding
facet_names <- list(
  'Up_inASD'="Increased TF binding", #"TF binding up in ASD"
  'Down_inASD'="Decreased TF binding" #"TF binding down in ASD"
)
facet_labeller <- function(variable,value){return(facet_names[value])}

TextSize = 5
PointSize = 1
LegendSize = 0.5

Fig2e = Npeaks_BoundByTF_EaSample_melt %>%
  ggplot(aes(x = TFmotif, y = Npeak_bound, col = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(), alpha = 0.5, size = PointSize/2) + 
  theme_bw() +
  ylab("TF-bound open chromatin regions") +  # "Number of ATAC peaks bound by TF"
  xlab("TF motifs") +
  ggforce::facet_row(vars(change), scales = 'free', space = 'free', labeller = facet_labeller) + # refer to https://stackoverflow.com/questions/52341385/how-to-automatically-adjust-the-width-of-each-facet-for-facet-wrap
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(size = TextSize),
        legend.text = element_text(size = TextSize), 
        title = element_text(size = TextSize),
        legend.margin = margin(0,unit = "lines"),
        legend.key.size = unit(LegendSize, 'lines'),
        strip.text.x = element_text(size = TextSize)
  )

pdf("../../plots/02_TF/2e_01_num_ATACpeaks_boundBy_TF_Fig2e.pdf", height = 2, width = 4.7) # original file similar to Figure2e_Number_of_ATACpeaks_BoundBy_TF_BindingScoreGoingInSameDirWithDGE.pdf
print(Fig2e)
dev.off()

# Observation:
# CTCF shows lower binding score and gene expression. All three CTCF motifs trend towards binding lower number of ATAC peaks bound in ASD. -> PlotAggregate shows that on either ASD bound or CTL bound regions, CTCF binding is lower in ASD vs. CTL.
