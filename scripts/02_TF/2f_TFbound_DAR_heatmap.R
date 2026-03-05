### Goal: This script focuses on # 3) heatmap of ATAC log(CPM+5) in ASD and CTL samples (row x col = peak x sample, left: all CTL samples, right: all ASD samples) 

rm(list = ls())
options(stringsAsFactors = F)

library(tidyverse)
library(ggplot2)
library(ggforce)
#install.packages("UpSetR")
library(UpSetR)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/") # ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/ind_BINDetect_Jaspar2022/

########## pipeline ########## 
### 1. Load ATAC sample covariates and TFs of interest
### 2. Catalog TF-bound ATAC peaks in any samples
### 3. Bootstrap on the number of DARs bound by each TFoi and save plots
### 4. Plot the histograms for each TFoi
##############################

#lnames = load("2d_diffTOBIAS_DGE.rda") # "DiffBS_lm_Volcano" "TFclusters". # original file: DiffBS_lm_Volcano_wCluster_DGE.rda
#lnames = load("2e_01_Npeaks_BoundByTF_EaSample.rda") # "Npeaks_BoundByTF_EaSample"     "Npeaks_BoundByTF_EaSample_melt". # original file: Npeaks_BoundByTF_EaSample.rda
#lnames = load("2e_01_num_ATACpeaks_boundBy_TFoi.rda") # Npeaks_BoundByTF_EaSample_melt. # original file: Number_of_ATACpeaks_BoundBy_TF_BindingScoreGoingInSameDirWithDGE.rda

## Upset showing the overlap of ATAC peaks bound by members of the FOS/JUN/BATF/NFE2 TF cluster

lnames = load("2d_TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda") # original file: TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda
TFoi = TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF; rm(TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF)

tf_bound = vector("list")
for (tf in c("BATF", "FOSL1", "BACH1", "NFE2")) {
  tfoi = TFoi[grepl(tf, TFoi)]
  bound = read.table(paste0("TFoi_bound_ATACpeaks/",tfoi, "_BoundInAnySample.txt"), header = T)
  bound$peak = paste0(bound$peak_chr, "_", bound$peak_start, "_", bound$peak_end)
  #assign(paste0(tf,"_bound"), unique(bound$peak))
  tf_bound[[tf]] = unique(bound$peak)
}
save(tf_bound, file = "UpSetR_OverlapOfNumATACBroundBy_Bach1BatfFosl1Nfe2_TFcluster.rda")

pdf("UpSetR_OverlapOfNumATACBroundBy_Bach1BatfFosl1Nfe2_TFcluster.pdf", width = 8, height = 5, onefile = F)
upset(fromList(tf_bound), order.by = "freq")
dev.off()
# Observation: 
# Over 20k ATAC peaks are co-bound by BACH1/FOSL1/BATF.

tf_bound2 = vector("list")
tfoi = TFoi[grepl("CTCF", TFoi)]
for (tfoi_cur in tfoi) {
  bound = read.table(paste0("TFoi_bound_ATACpeaks/",tfoi_cur, "_BoundInAnySample.txt"), header = T)
  bound$peak = paste0(bound$peak_chr, "_", bound$peak_start, "_", bound$peak_end)
  #assign(paste0(tf,"_bound"), unique(bound$peak))
  tf_bound2[[tfoi_cur]] = unique(bound$peak)
}
save(tf_bound2, file = "UpSetR_OverlapOfNumATACBroundBy_CTCF_TFcluster.rda")

pdf("UpSetR_OverlapOfNumATACBroundBy_CTCF_TFcluster.pdf", width = 8, height = 5, onefile = F)
upset(fromList(tf_bound2), order.by = "freq")
dev.off()
# Observation: 
# Around half of ATAC peaks bound by each CTCF motif (~22k) are co-bound by all 3 CTCF motifs.

# So for heatmaps, I only need to show change of a. CTCF all 3 motifs co-bound b.co-bound by FOSL1, BATF, BACH1 ATAC peaks.

## heatmap
lnames = load("UpSetR_OverlapOfNumATACBroundBy_Bach1BatfFosl1Nfe2_TFcluster.rda")
lnames = load("UpSetR_OverlapOfNumATACBroundBy_CTCF_TFcluster.rda")
lnames = load("../../5_DiffATAC/CQN.rda") # "RPKM.cqn"   "cqn.fit"    "Covariates" "peak_info"  "counts_tbl"
# original
#lnames = load("../../12_ARI/ATACregressed_ExceptDiagnosisBAregionAge_pARI.rda") # "ATAC_regressed" "p_ARI" "fdr_ARI"
# re-run to make figures for manuscript
lnames = load("../../12_ARI/v14_ATACregressed_ExceptDiagnosisBAregionAgeBatch.rda") # "ATAC_regressed"
# v2_ATACregressed_ExceptDiagnosisBAregionAge.rda doesn't give better heatmap.

CTCF_bound1 = read.table(paste0("TFoi_bound_ATACpeaks/CTCF-MA1930.1_BoundInAnySample.txt"), header = T)
CTCF_bound2 = read.table(paste0("TFoi_bound_ATACpeaks/CTCF-MA0139.1_BoundInAnySample.txt"), header = T)
CTCF_bound3 = read.table(paste0("TFoi_bound_ATACpeaks/CTCF-MA1929.1_BoundInAnySample.txt"), header = T)
tmp = inner_join(CTCF_bound1, CTCF_bound2)
CTCF_bound = inner_join(tmp, CTCF_bound3) # 22242, match the number on UpSetR
colnames(CTCF_bound) = c("chr", "start", "end")
CTCF_bound = left_join(CTCF_bound, DiffATAC)
CTCF_bound_ATACdown = CTCF_bound %>%
  filter(ATAC_logFC < 0, ATAC_FDR < 0.05) # 523
CTCF_bound_ATACdown$ATACpeak = paste0(CTCF_bound_ATACdown$chr, "_", CTCF_bound_ATACdown$start, "_", CTCF_bound_ATACdown$end)

BACH1_bound = read.table(paste0("TFoi_bound_ATACpeaks/BACH1-MA1633.2_BoundInAnySample.txt"), header = T)
FOSL1_bound = read.table(paste0("TFoi_bound_ATACpeaks/FOSL1-MA0477.2_BoundInAnySample.txt"), header = T)
BATF_bound = read.table(paste0("TFoi_bound_ATACpeaks/BATF-MA1634.1_BoundInAnySample.txt"), header = T)
tmp = inner_join(BACH1_bound, FOSL1_bound)
co_bound = inner_join(tmp, BATF_bound) # 21414, match the number on UpSetR
colnames(co_bound) = c("chr", "start", "end")
co_bound = left_join(co_bound, DiffATAC)
co_bound_ATACup = co_bound %>%
  filter(ATAC_logFC > 0, ATAC_FDR < 0.05) # 670
co_bound_ATACup$ATACpeak = paste0(co_bound_ATACup$chr, "_", co_bound_ATACup$start, "_", co_bound_ATACup$end)

# CTCF_bound_ATACdown_RPKMcqn = RPKM.cqn[rownames(RPKM.cqn) %in% CTCF_bound_ATACdown$ATACpeak,]
# co_bound_ATACup_RPKMcqn = as.data.frame(RPKM.cqn[rownames(RPKM.cqn) %in% co_bound_ATACup$ATACpeak,])
CTCF_bound_ATACdown_RPKMcqn = ATAC_regressed[rownames(ATAC_regressed) %in% CTCF_bound_ATACdown$ATACpeak,]
co_bound_ATACup_RPKMcqn = as.data.frame(ATAC_regressed[rownames(ATAC_regressed) %in% co_bound_ATACup$ATACpeak,])

Covariates$Diagnosis = factor(Covariates$Diagnosis, levels = c("CTL", "ASD"))
Covariates$Cortex = factor(Covariates$Cortex, levels = c("Parietal", "Frontal", "Temporal"))
Covariates = Covariates %>%
  arrange(Diagnosis, Cortex, Age)

CTCF_bound_ATACdown_RPKMcqn = CTCF_bound_ATACdown_RPKMcqn[,Covariates$sample_name]
co_bound_ATACup_RPKMcqn = co_bound_ATACup_RPKMcqn[,Covariates$sample_name]

# install.packages("pheatmap")
library(pheatmap)

annoBar = Covariates[,c("Diagnosis", "Age", "Cortex")] # , "Batch"
annoCol = list(Diagnosis = c(CTL = "aquamarine", ASD = "#F8766D"), 
               Cortex = c(Frontal = "dodgerblue", Parietal = "gold", Temporal = "plum2"), # Frontal = "#00BA38"
               Sex = c(M = "#C77CFF", F = "pink"))

# hm_CTCF = pheatmap::pheatmap(CTCF_bound_ATACdown_RPKMcqn, scale = "row", show_colnames = F, show_rownames = F, annotation_col = annoBar, annotation_colors = annoCol, main = "Normalized ATAC-seq log2(CPM)") # scale = "row"
hm_CTCF = pheatmap::pheatmap(CTCF_bound_ATACdown_RPKMcqn, scale = "row", show_colnames = F, show_rownames = F, cluster_cols = F, annotation_col = annoBar, annotation_colors = annoCol, main = "Normalized ATAC-seq log2(CPM)", treeheight_row = 0) # scale = "row"

hm_FOSBACH1BATF = pheatmap::pheatmap(co_bound_ATACup_RPKMcqn, scale = "row", show_colnames = F, show_rownames = F, cluster_cols = F, annotation_col = annoBar, annotation_colors = annoCol, main = "Normalized ATAC-seq log2(CPM)", treeheight_row = 0) # scale = "row"

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
require(RColorBrewer)
require(circlize)
require(digest)
require(cluster)

myCol <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
col_fun = colorRamp2(c(5, 33, 60), brewer.pal(n = 3, name = "Greens"))
annoCol2 = list(Diagnosis = c(CTL = "aquamarine", ASD = "#F8766D"), 
               Cortex = c(Frontal = "dodgerblue", Parietal = "gold", Temporal = "plum2"), # Frontal = "#00BA38"
               Sex = c(M = "#C77CFF", F = "pink"),
               Age = col_fun)
hmAnno_show = HeatmapAnnotation(
  df = annoBar, which = 'col', # 'col' (samples) or 'row' (gene) annotation
  na_col = "white", col = annoCol2, show_annotation_name = T) # , simple_anno_size = unit(2, "mm"), gp = gpar(col = "grey", lwd = 0.5))

scaled_CTCF = t(scale(t(CTCF_bound_ATACdown_RPKMcqn)))
hm_CTCF2 = Heatmap(scaled_CTCF, col = myCol, name = "Normalized ATAC-seq log2(CPM)", top_annotation = hmAnno_show, heatmap_legend_param = list(legend_direction = "horizontal"), cluster_columns = F, show_row_dend = F , show_row_names = F, show_column_names = F, column_title = paste0("CTCF-bound ATAC Peaks down-regulated in ASD (N=", formatC(nrow(CTCF_bound_ATACdown_RPKMcqn), format = "d", big.mark = ","), ")"), column_title_gp = gpar(fontsize = 12)) #, row_split = 5, column_split = 5, row_title_gp = gpar(fontsize = 8)) # row_dend_width = unit(0.5, "cm"), 

n_CTL = length(which(Covariates$Diagnosis == "CTL"))

pdf("Heatmap_CTCF_bound_ASDdownATAC.pdf", width = 8, height = 5)
draw(hm_CTCF2, heatmap_legend_side="bottom", annotation_legend_side="right", legend_grouping = "original") # , newpage = F
#draw(hm_CTCF2, ht_gap = unit(3, "mm"))
decorate_heatmap_body("Normalized ATAC-seq log2(CPM)", {
  grid.lines(c(n_CTL/nrow(Covariates),n_CTL/nrow(Covariates)), c(0,1), gp = gpar(lwd = 2,lty = 2, col = "black"))
  })
dev.off()

save(hm_CTCF2, n_CTL, Covariates, file = "Heatmap_CTCF_bound_ASDdownATAC.rda")

scaled_FOSBACH1BATF = t(scale(t(co_bound_ATACup_RPKMcqn)))
hm_FOSBACH1BATF2 = Heatmap(scaled_FOSBACH1BATF, col = myCol, name = "Normalized ATAC-seq log2(CPM)", top_annotation = hmAnno_show, heatmap_legend_param = list(legend_direction = "horizontal"), cluster_columns = F, show_row_dend = F , show_row_names = F, show_column_names = F, column_title = paste0("FOSL1/BACH1/BATF co-bound ATAC Peaks up-regulated in ASD (N=", formatC(nrow(co_bound_ATACup_RPKMcqn), format = "d", big.mark = ","), ")"), column_title_gp = gpar(fontsize = 12)) #, row_split = 5, column_split = 5, row_title_gp = gpar(fontsize = 8)) # row_dend_width = unit(0.5, "cm"), 

pdf("Heatmap_FOSBACH1BATF_cobound_ASDupATAC.pdf", width = 8, height = 5)
draw(hm_FOSBACH1BATF2, heatmap_legend_side="bottom", annotation_legend_side="right", legend_grouping = "original") # , newpage = F
#draw(hm_CTCF2, ht_gap = unit(3, "mm"))
decorate_heatmap_body("Normalized ATAC-seq log2(CPM)", {
  grid.lines(c(n_CTL/nrow(Covariates),n_CTL/nrow(Covariates)), c(0,1), gp = gpar(lwd = 2,lty = 2, col = "black"))
})
dev.off()

save(hm_FOSBACH1BATF2, n_CTL, Covariates, file = "Heatmap_FOSBACH1BATF_cobound_ASDupATAC.rda")
# Is there evidence that regional differential ATAC is attenuated in ASD?

## Next step: link ATAC peaks to genes, show change in gene expression. Show change in Hi-C loops. 
# And all the results section in Corces 2018 Science

