### Goal: This script focuses on # 2) Randomly draw n ATAC peaks -> build null distribution of probability of ATAC diff in the same dir -> p-val of observed probablity of ATAC bound being diff in the same dir 

rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(reshape2)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/") # original directory: /u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/ind_BINDetect_Jaspar2022/all_bindetect_results

########## pipeline ########## 
### 1. Load TF binding and DAR data
### 2. Catalog TF-bound ATAC peaks in any samples
### 3. Bootstrap on the number of DARs bound by each TFoi and save plots
### 4. Plot the histograms for each TFoi
##############################

### 1. Load TF binding and DAR data
lnames = load("2d_TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda") # original file: TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda
TFoi = TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF; rm(TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF)
lnames = load("2d_diffTOBIAS_DGE.rda") # "DiffBS_lm_Volcano" "TFclusters". # original file: DiffBS_lm_Volcano_wCluster_DGE.rda
lnames = load("../01_DAR/1a_03_DiffATAC.rda") # DiffATAC. original file: ../../5_DiffATAC/DiffATAC.rda 

### 2. Catalog TF-bound ATAC peaks in any samples
# This part 2 is modified from original script 10_3_7_Link_DiffTFbsAndge_to_ATACchange.R
DIR = paste0(getwd(),"2c_01_TOBIAS/ind_BINDetect/") # original directory: "/u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/10_TOBIAS/ind_BINDetect/all_bindetect_results" 

for (tf in TFoi) {
  print(tf)
  peak_bound = matrix(nrow = 0, ncol = 3)
  for (s in Covariates$sample_name) {
    TF_fullname = rownames(DiffBS_lm_Volcano)[DiffBS_lm_Volcano$TFmotif == tf]
    TF_dir = paste0(TF_fullname,"_",TF_fullname)
    bedfile = paste0(DIR, s,"/",TF_dir,"/",TF_dir,"_overview.txt") # original start with "../"
    bed = read.table(bedfile, header = T)
    bed = bed[bed[,ncol(bed)] == 1,]
    peak_bound_s = unique(bed[,7:9])
    peak_bound = unique(rbind(peak_bound, peak_bound_s))
  }
  print(dim(peak_bound))
  write.table(peak_bound, file = paste0("2e_02_TFoi_bound_ATACpeaks/", tf, "_BoundInAnySample.txt"), quote = F, row.names = F, col.names = T) # originally no TFoi_bound_ATACpeaks/
}
# Prints:
# E2F1-MA0024.2      20912
# BATF_HUMAN.H10MO.S 31155
# FOSL1-MA0477.1     26271 
# JUN_HUMAN.H10MO.A  31306
# NFE2-MA0841.1      19559
# NFIL3-MA0025.1     13111
# CTCF-MA0139.1      49433

### 3. Bootstrap on the number of DARs bound by each TFoi and save plots
TextSize = 5
PointSize = 1
LegendSize = 0.5
Lwd = 0.2

plot_diffATAC = list()

for (tf in TFoi) {
  print(tf)
  
  ## Number of diffATAC in the same dir among the TF-bound ATAC peaks
  bound = read.table(paste0(DIR, "2e_02_TFoi_bound_ATACpeaks/",tf, "_BoundInAnySample.txt"), header = T) # original directory: TFoi_bound_ATACpeaks/
  colnames(bound) = c("chr", "start", "end")
  bound = left_join(bound, DiffATAC)
  Obs_up = length(which(bound$ATAC_logFC > 0 & bound$ATAC_FDR < 0.05))
  Obs_down = length(which(bound$ATAC_logFC < 0 & bound$ATAC_FDR < 0.05))
  Odds_up = (Obs_up/nrow(bound))/(1-Obs_up/nrow(bound))
  Odds_down = (Obs_down/nrow(bound))/(1-Obs_down/nrow(bound))
  
  ## N ATAC peaks to draw
  Ntodraw = nrow(bound)
  Type = TFclusters$change[TFclusters$TFmotif == tf] # Down_inASD
  
  set.seed(183)
  Ndiff_up = c(); OR_up_bootstrap = c()
  Ndiff_down = c(); OR_down_bootstrap = c()
  for (B in 1:1000) {
    idx = sample(1:nrow(DiffATAC), size = Ntodraw)
    idx_ATAC = DiffATAC[idx,]
    n_diff_up = length(which(DiffATAC$ATAC_logFC[idx] > 0 & DiffATAC$ATAC_FDR[idx] < 0.05))
    n_diff_down = length(which(DiffATAC$ATAC_logFC[idx] < 0 & DiffATAC$ATAC_FDR[idx] < 0.05))
    Ndiff_up = c(Ndiff_up, n_diff_up)
    Ndiff_down = c(Ndiff_down, n_diff_down)
    
    Odds_up_current = (n_diff_up/nrow(bound))/(1-n_diff_up/nrow(bound))
    Odds_down_current = (n_diff_down/nrow(bound))/(1-n_diff_down/nrow(bound))
    OR_up_current = Odds_up/Odds_up_current
    OR_down_current = Odds_down/Odds_down_current
    OR_up_bootstrap = c(OR_up_bootstrap, OR_up_current)
    OR_down_bootstrap = c(OR_down_bootstrap, OR_down_current)
  }
  
  ## histogram, p-val, OR
  p_up = length(which(Ndiff_up > Obs_up))/length(Ndiff_up)
  p_down = length(which(Ndiff_down > Obs_down))/length(Ndiff_down)
  OR_up = mean(OR_up_bootstrap)
  OR_down = mean(OR_down_bootstrap)
  print(paste0(tf, ", p_up = ", p_up, ", p_down = ", p_down, 
               ", OR_up = ", OR_up, ", OR_down = ", OR_down))
  
  df = as.data.frame(cbind(Ndiff_up, Ndiff_down))
  
  p_cur1 = df %>%
    ggplot(aes(x = Ndiff_up)) +
    geom_histogram(col = "black", fill = "#F8766D", alpha = 0.2) + # red
    geom_vline(xintercept = Obs_up, col = "red", size = Lwd*5) +
    theme_bw() +
    xlab("Regulatory regions with increased accessibility") +
    ylab("Frequency") +
    ggtitle(paste0(tf, " (", Ntodraw," regions)\nOR = ", round(OR_up,2),  ifelse(p_up > 0, paste0(", p = ", round(p_up,2)), ", p < 1e-3"), " (one-sided)")) +
    theme(
      text = element_text(size = TextSize),
      axis.text = element_text(size = TextSize),
      legend.text = element_text(size = TextSize), # required to have legends at 5 pt
      title = element_text(size = TextSize),
      plot.title = element_text(hjust = 0.5, size = TextSize)
    )
  
  p_cur2 = df %>%
    ggplot(aes(x = Ndiff_down)) +
    geom_histogram(col = "black", fill = "#00BFC4", alpha = 0.2) + # blue
    geom_vline(xintercept = Obs_down, col = "red", size = Lwd*5) +
    theme_bw() +
    xlab("Regulatory regions with decreased accessibility") +
    ylab("Frequency") +
    ggtitle(paste0(tf, " (", Ntodraw," regions)\nOR = ", round(OR_down,2),  ifelse(p_down > 0, paste0(", p = ", round(p_down,2)), ", p < 1e-3"), " (one-sided)")) +
    theme(
      text = element_text(size = TextSize),
      axis.text = element_text(size = TextSize),
      legend.text = element_text(size = TextSize),
      title = element_text(size = TextSize),
      plot.title = element_text(hjust = 0.5, size = TextSize))
  
  plot_diffATAC[[paste0(tf, "_Up")]] = p_cur1
  plot_diffATAC[[paste0(tf, "_Down")]] = p_cur2
}
names(plot_diffATAC)
save(plot_diffATAC, file = "2e_02_Histograms_num_DARs_boundByTF.rda") # original file: FigureSup5e_Histograms_TFbound_NdiffATAC.rda

### 4. Plot the histograms for each TFoi
# The following scripts are adapted from Rscripts_v5_DanEdit/SupFigure5.R
rm(list = ls())
library(ggpubr)
library(purrr)

TextSize = 5
PointSize = 1
LegendSize = 0.5
Lwd = 0.2

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/") 
lnames = load("2e_02_Histograms_num_DARs_boundByTF.rda")

pdf("../../plots/02_TF/2e_02_Histograms_num_DARs_boundByTF_FigSup5e_Fig2f.pdf", width = 7, height = 7) # originally Figures_v5/Supplementary Figure 5e_Histograms_TFbound_NdiffATAC.pdf in 
FigSup5e = ggarrange(plotlist = rev(plot_diffATAC), nrow = 4, ncol = 4)
print(FigSup5e) # originally Fig_Sup4b
dev.off()

## To get Main Figure 2f:
library(cowplot)
Fig2f = plot_grid(plot_diffATAC[["CTCF-MA0139.1_Down"]], plot_diffATAC[["BATF-MA1634.1_Up"]], plot_diffATAC[["NFIL3-MA0025.2_Up"]], nrow = 1, ncol = 3)

