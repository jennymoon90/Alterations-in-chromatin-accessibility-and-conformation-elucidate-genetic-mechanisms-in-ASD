rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")

########## pipeline ########## 
### 1. Load CHART-assigned DAR target genes
### 2. Quantify the number of up/down-reg DARs at promoter and distal regions
### 3. Plot bargraph
##############################

### 1. Load CHART-assigned DAR target genes
lnames = load("1c_01_CHART_assigned_distalATAC_gene_pairs.rda") # "sig_CorATAC"      "sig_CorATAC_uniq" "promATAC"
# original file: ../v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda

DARprom_uniq = unique(promATAC[which(promATAC$ATAC_FDR < 0.05), c("ATACpeak", "ATAC_logFC", "ATAC_FDR")])
(Prom_total = nrow(DARprom_uniq)) # 1305
(Prom_up = length(which(DARprom_uniq$ATAC_logFC > 0))) # 799 ~ 61%
(Prom_down = length(which(DARprom_uniq$ATAC_logFC < 0))) # 506 ~ 39%

DARdist_uniq = unique(sig_CorATAC_uniq[which(sig_CorATAC_uniq$distalATAC_FDR < 0.05), c("distal_ATAC", "distalATAC_logFC", "distalATAC_FDR")])
(DistalA_total = nrow(DARdist_uniq)) # 1043
(DistalA_up = length(which(DARdist_uniq$distalATAC_logFC > 0))) # 420 ~ 40%
(DistalA_down = length(which(DARdist_uniq$distalATAC_logFC < 0))) # 623 ~ 60%

### 2. Quantify the number of up/down-reg DARs at promoter and distal regions
counts = data.frame(
  ATAC_type = rep(c("Promoter", "Distal_assigned"), each = 2), 
  dysreg = rep(c("Down-reg", "Up-reg"), 2), 
  Freq = c(Prom_down, Prom_up, DistalA_down, DistalA_up),
  Percent = c(Prom_down/Prom_total, Prom_up/Prom_total, DistalA_down/DistalA_total, DistalA_up/DistalA_total))

counts$dysreg = factor(counts$dysreg, levels = c("Up-reg", "Down-reg"))
counts$ATAC_type = factor(counts$ATAC_type, levels = c("Promoter", "Distal_assigned"))

### 3. Plot bargraph
TextSize = 5 
AnnoSize = 1.8
LegendSize = 0.5

Fig1g = counts %>% 
  ggplot(aes(x = ATAC_type, y = Percent, fill = dysreg)) +
  geom_col(position = position_stack()) +
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = AnnoSize) +   theme_bw() +
  labs(fill = "DAR direction\nin ASD") + # 
  xlab("Location of DARs relative to target genes") +
  scale_x_discrete(labels=c('Promoter', 'Distal')) + # , 'Distal\nassigned', 'Distal not\nassigned'
  ggtitle("Number of differentially accessible regions (DARs)") +
  theme(text = element_text(size = TextSize),
        axis.text = element_text(size = TextSize),
        legend.text = element_text(size = TextSize), # required to have legends at 5 pt
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.4, size = TextSize),
        legend.key.size = unit(LegendSize, 'lines'),
        legend.margin = margin(0,unit = "lines"))

pdf("../../plots/01_DAR/1c_02_num_DARs_PromoterVsDistal_Fig1g.pdf", width = 2, height = 2) 
# original file: ~/Documents/Documents/Geschwind_lab/LAB/Manuscripts/ASD_Dup15q_project/Figures_v3/Figure1g_diff_ATAC_Type_PromoterDistalAssignedOrNot_PromTSSminus2kbplus100bp_DistalcorWithPromTSSminusplus2kb.pdf
print(Fig1g)
dev.off()

