rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(reshape2)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/") # original: ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/10_TOBIAS/ind_BINDetect_Jaspar2022/

########## pipeline ########## 
### 1. Load TF binding scores by sample
### 2. Calculate the change in TF binding score (ASD vs CTL)
### 3. Cluster by TF family
### 4. Add TF gene expression
### 5. Scatter plot of TF binding score and the TF gene expression
### 6. Save the significant TFs
##############################

### 1. Load TF binding scores by sample
lnames = load("2c_02_Jaspar2022_TFmeanScore_by_sample.rda") # df. # original: Jaspar2022hs_TFmotifs_mean_score_by_sample.rda
colnames(df) = sub("_mean.score", "", colnames(df))
colnames(df) = sub("\\.", "-", colnames(df))
colnames(df)[2:5] = sub("X", "", colnames(df)[2:5])

## Filter for the motifs of interest
lnames = load("2c_03_diffTOBIAS_volcano.rda") # DiffBS_lm_Volcano. # original file: DiffBS_lm_Volcano.rda
df = df[df$motif_id %in% rownames(DiffBS_lm_Volcano)[DiffBS_lm_Volcano$TFlabel != ""],] # 45 TFmotifs
df$TFmotif = DiffBS_lm_Volcano$TFmotif[match(df$motif_id, rownames(DiffBS_lm_Volcano))]
df = df[,-1]

### 2. Calculate the change in TF binding score (ASD vs CTL)
## Melt motif binding score by sample dataframe (df -> df_melt)
df_melt = melt(df, variable.name = "sample_name", value.name = "mean_binding_score")

## Attach sample meta data (Age, Diagnosis)
lnames = load("../../data/ATAC_library_covariates.rda") # df_wAllCov. # original: ../../3_QC/BiolTechCov_and_QC.rda
colnames(df_wAllCov)
datMeta = df_wAllCov[match(df_melt$sample_name, df_wAllCov$sample_name), c("sample_name", "Diagnosis", "Age", "Sex", "Cortex", "Region")]
df_melt = cbind(df_melt, datMeta[,-1])

## Calculate the change
df_melt$change = DiffBS_lm_Volcano$col[match(df_melt$TFmotif, DiffBS_lm_Volcano$TFmotif)]

### 3. Cluster by TF family
df_melt$cluster = DiffBS_lm_Volcano$cluster[match(df_melt$TFmotif, DiffBS_lm_Volcano$TFmotif)]

TFclusters = unique(df_melt[,c("TFmotif", "cluster", "change")])
TFclusters = TFclusters %>%
  arrange(change, cluster)
TFclusters$TF_shortname = DiffBS_lm_Volcano$TF_shortname[match(TFclusters$TFmotif, DiffBS_lm_Volcano$TFmotif)]

TFclusters$TFs = TFclusters$TF_shortname
for (i in 1:nrow(TFclusters)) {
  c = TFclusters$cluster[i]
  idx = which(TFclusters$cluster == c)
  if (length(idx) > 1) {
    TFclusters$TFs[i] = paste0(unique(TFclusters$TF_shortname[idx]), collapse = "/")
  }
}

df_melt$TFs = TFclusters$TFs[match(df_melt$TFmotif, TFclusters$TFmotif)]
unique(df_melt$TFs)
df_melt$TFs[which(df_melt$TFs == "BACH1/BACH2/BATF/BATF::JUN/BATF3/BNC2/FOS/FOS::JUN/FOS::JUNB/FOS::JUND/FOSB::JUNB/FOSL1/FOSL1::JUN/FOSL1::JUNB/FOSL1::JUND/FOSL2/FOSL2::JUN/FOSL2::JUNB/FOSL2::JUND/JDP2/JUN::JUNB/JUNB/JUND/NFE2")] = "BACH(1/2)/BATF(3)/BNC2/FOS(L1/2)/JUN(B/D)/JDP2/NFE2"
df_melt$TFs[which(df_melt$TFs == "MEF2A/MEF2B/MEF2C/MEF2D")] = "MEF2A/B/C/D"
df_melt$TFs[which(df_melt$TFs == "BHLHA15/BHLHE23/OLIG1/OLIG2/OLIG3")] = "BHLHA15/BHLHE23/OLIG(1/2/3)"

idx_dup = which(duplicated(TFclusters[,c("cluster", "TFs")]))
TF_dups = TFclusters$TFmotif[idx_dup] # 33
df_melt2 = df_melt[-which(df_melt$TFmotif %in% TF_dups),]
df_melt2 = df_melt2 %>%
  arrange(change, TFs)
df_melt2$TFs = factor(df_melt2$TFs, levels = unique(df_melt2$TFs))

### 4. Add TF gene expression
library(readxl)
Jill_DEG = read_excel("../../data/DEG_Gandal2022Nature.xlsx") # 24836 genes in DEG data from Jill's pan-cortical study. # original file: ~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx

DiffBS_lm_Volcano$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(DiffBS_lm_Volcano$TF_shortname, Jill_DEG$external_gene_name)]
DiffBS_lm_Volcano$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(DiffBS_lm_Volcano$TF_shortname, Jill_DEG$external_gene_name)]

TFclusters = left_join(TFclusters, DiffBS_lm_Volcano[,c("TFmotif", "magnitude", "p", "sig", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR")])

TFclusters$GE_sig = ifelse(TFclusters$WholeCortex_ASD_FDR < 0.05, "significant", "n.s.") # use FDR < 0.05 rather than 0.1! If 0.1, one more point in Q2.
TFclusters$GE_sig[is.na(TFclusters$GE_sig)] = "n.s."
TFclusters$GE_sig = factor(TFclusters$GE_sig, levels = c("significant", "n.s."))

TFclusters$GE_label = ifelse(TFclusters$GE_sig == "significant" & TFclusters$magnitude * TFclusters$WholeCortex_ASD_logFC > 0, TFclusters$TF_shortname, "")
TFclusters$GE_label = ifelse(TFclusters$magnitude * TFclusters$WholeCortex_ASD_logFC > 0, TFclusters$TF_shortname, "")

save(DiffBS_lm_Volcano, TFclusters, file = "2c_05_diffTOBIAS_DGE.rda") # original file: DiffBS_lm_Volcano_wCluster_DGE.rda

### 5. Scatter plot of TF binding score and the TF gene expression
TextSize = 5
AnnoSize = 1.8
PointSize = 1
LegendSize = 0.5
Lwd = 0.2

Fig2d = TFclusters %>%
  ggplot(aes(x = WholeCortex_ASD_logFC, y = magnitude, col = GE_sig, label = GE_label)) +
  geom_point(alpha = 0.7, size = PointSize/2) +
  geom_text_repel(box.padding = Lwd, max.overlaps = Inf, size = AnnoSize, seed = 42, segment.size = Lwd, segment.alpha = 0.5) +
  ylim(-0.013, 0.015) +
  xlim(-0.25, 1) +
  labs(col = "Differential TF gene expression") +
  theme_bw() +
  xlab("Differential TF gene expression\n(logFC, ASD vs CTL)") + 
  ylab("TF binding score (ASD - CTL)") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(text = element_text(size = TextSize),
        axis.text = element_text(size = TextSize),
        legend.text = element_text(size = TextSize), 
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.5, size = TextSize),
        legend.margin = margin(0,0,-0.3,0,unit = "lines"),
        legend.key.size = unit(LegendSize, 'lines'),
        legend.position = "top"
  )

pdf("../../plots/02_TF/2c_05_scatter_diffTOBIAS_DGE_Fig2d.pdf", width = 2.5, height = 2)
print(Fig2d)
dev.off()

### 6. Save the significant TFs
(TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF = TFclusters$TFmotif[which(TFclusters$GE_sig == "significant" & TFclusters$magnitude * TFclusters$WholeCortex_ASD_logFC > 0)]) # "BACH1"  "BATF"  "FOSL1" "NFE2"  "NFIL3" "xx::xx"
TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF = c(TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF, TFclusters$TFmotif[TFclusters$TF_shortname == "CTCF"])

save(TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF, file = "2c_05_TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda") # original file: TFmotif_BindingScore_and_DGE_BothSig_InSameDir_plusCTCF.rda


