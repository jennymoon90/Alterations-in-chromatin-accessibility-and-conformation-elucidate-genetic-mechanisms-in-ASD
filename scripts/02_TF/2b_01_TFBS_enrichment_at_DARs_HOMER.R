rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(Repitools)
library(reshape2)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/")

########## pipeline ########## 
### 1. Prepare bed files for HOMER TF enrichment
### 2. Run HOMER
### 3. Organize TF binding motif enrichment results
### 4. Barplot
##############################

### 1. Prepare bed files for HOMER TF enrichment

#setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/")

## Load DARs
lnames = load("../01_DAR/1a_03_DiffATAC.rda") # DiffATAC. # original file: ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/DiffATAC.rda

## Filter to differential ATAC peaks
ATAC_up = DiffATAC[which(DiffATAC$ATAC_logFC > 0 & DiffATAC$ATAC_FDR < 0.05),]
ATAC_down = DiffATAC[which(DiffATAC$ATAC_logFC < 0 & DiffATAC$ATAC_FDR < 0.05),]

## Write to bed files
# original file in ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/
write.table(ATAC_up, file = paste0("2a_01_ATACup.bed"), quote = F, row.names = F, col.names = F, sep = "\t") # original file: upATAC_fdrthreshold0.05logfcthreshold0.bed
write.table(ATAC_down, file = paste0("2a_01_ATACdown.bed"), quote = F, row.names = F, col.names = F, sep = "\t") # original file: downATAC_fdrthreshold0.05logfcthreshold0.bed

### 2. Run HOMER: 2a_02_HOMER_MotifScan_in_DARupORdown_Jaspar2022.sh
# original script: 9_4_HOMER_MotifScan_in_DiffATACupORdown_Jaspar2022.sh

# For motif, used jaspar2022_BrExpTF.motifs generated using 14_10_TFBS_enrichment_in_4dif_categories_of_gene_BasedOn_Nloop_TPM_Jaspar2022.R and 14_11_HOMER_scan_motif_in_promoterATAC_4categories_Jaspar2022.sh (this is the JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt -> restrict to human brain expressed TFs)

### 3. Organize TF binding motif enrichment results
rm(list = ls())

(motif_files = list.files(pattern = "*_knownResults.txt")) # 2 files. original directory: ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/Jaspar2022_results/

for (motif_file in motif_files) {
  (new_name = substr(motif_file,1,nchar(motif_file) - 4))
  new_name = gsub("2a_01_", "", new_name)
  new_name = gsub("_Jaspar2022", "", new_name)
  
  new_motif = read.table(paste0(motif_file))
  colnames(new_motif) = c("Consensus", "P_value", "logP_value", "q_value_Benjamini", 
                          "Number_of_Target_Sequences_with_Motif","Percent_of_Target_Sequences_with_Motif",
                          "Number_of_Background_Sequences_with_Motif","Percent_of_Background_Sequences_with_Motif")
  new_motif$Motif_name = rownames(new_motif)
  assign(new_name, new_motif)
}

## Generate a matrix comparing the q-value of each motif.
ref = ATACdown_knownResults
ATACup_knownResults = ATACup_knownResults[match(rownames(ref), rownames(ATACup_knownResults)),]

Motif_comparison_main = 
  data_frame(Motif_name = ref$Motif_name,
             Consensus = ref$Consensus,
             
             pval_down = ATACdown_knownResults$P_value,
             pval_up = ATACup_knownResults$P_value,

             qval_down = ATACdown_knownResults$q_value_Benjamini,
             qval_up = ATACup_knownResults$q_value_Benjamini,
  )

idx = which(apply(Motif_comparison_main[,which(startsWith(colnames(Motif_comparison_main), "qval"))], 1, function(x) any(x < 0.05)))
Motif_comparison_main = 
  Motif_comparison_main[idx,] # only FDR<0.05 are shown

Motif_comparison_main$minus_log_qval_down = -log10(Motif_comparison_main$qval_down)
Motif_comparison_main$minus_log_qval_up = -log10(Motif_comparison_main$qval_up)

Motif_comparison_main$Motif_abbr = sapply(Motif_comparison_main$Motif_name, function(x) substr(x, 1, which(strsplit(x, "")[[1]] == "/")[1] -1))

Motif_comparison_main_melt = melt(Motif_comparison_main[,which(startsWith(colnames(Motif_comparison_main), "minus_log_qval") | colnames(Motif_comparison_main) == "Motif_abbr")])
colnames(Motif_comparison_main_melt)[3] = c("significance")
Motif_comparison_main_melt$variable = as.character(Motif_comparison_main_melt$variable)
Motif_comparison_main_melt$variable = gsub("minus_log_qval_", "", Motif_comparison_main_melt$variable)

Motif_comparison_main$sig_n = apply(Motif_comparison_main[,which(grepl("^qval", colnames(Motif_comparison_main)))], 1, function(x) length(which(x < 0.05)))
table(Motif_comparison_main$sig_n) # 3 motifs are enriched in both up and down-reg ATAC peaks, 35 motifs are uniquely enriched in either up or down-reg ATAC peaks.

Motif_comparison_main = Motif_comparison_main %>%
  arrange(dplyr::desc(sig_n), dplyr::desc(minus_log_qval_down), minus_log_qval_up, Motif_name)

duplicated(Motif_comparison_main$Motif_abbr) # BHLHE22 is duplicated as it has two motifs, but one is enriched in up-reg ATAC, while the other in down-reg ATAC. -> delete BHLHE22
Motif_comparison_main = Motif_comparison_main[Motif_comparison_main$Motif_abbr != "BHLHE22",]

Motif_comparison_main_melt$Motif_abbr =
  factor(Motif_comparison_main_melt$Motif_abbr,
         levels = Motif_comparison_main$Motif_abbr) 
Motif_comparison_main_melt = Motif_comparison_main_melt[!is.na(Motif_comparison_main_melt$Motif_abbr),]

Motif_comparison_main_melt$variable = gsub("down", "down-reg in ASD", Motif_comparison_main_melt$variable)
Motif_comparison_main_melt$variable = gsub("up", "up-reg in ASD", Motif_comparison_main_melt$variable)

## Save the results
save(Motif_comparison_main_melt, file = "2a_01_Motif_Enrichment_in_DARupORdown_HOMERq005_Jaspar2022.rda") # original file: /Volumes/DataTransferBwMac/Working_Dir/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/9_HOMER_diffATAC/Jaspar2022_results/Motif_Enrichment_in_DiffATACupORdown_HOMERq005_Jaspar2022.rda

### 4. Barplot
rm(list = ls())
lnames = load("2a_01_Motif_Enrichment_in_DARupORdown_HOMERq005_Jaspar2022.rda")

# Order TFBS by enrichment in Up-reg ATAC then Down-reg ATAC
Motif_comparison_main_melt$significance_Inf5 = ifelse(is.infinite(Motif_comparison_main_melt$significance), 5, Motif_comparison_main_melt$significance)
df = dcast(Motif_comparison_main_melt[,-3], Motif_abbr ~ variable)
df$delta = df$`up-reg in ASD` - df$`down-reg in ASD`
df$Motif_abbr = as.character(df$Motif_abbr)
TF_order = unique(df$Motif_abbr)
Motif_comparison_main_melt$Motif_abbr = factor(Motif_comparison_main_melt$Motif_abbr, levels = TF_order)

TextSize = 5
LegendSize = 0.5

Fig2b = Motif_comparison_main_melt %>% # originally Fig3a in Rscripts_v11/Figure2.R
  ggplot(aes(x = Motif_abbr, y = significance, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Enrichment of TFBS in differential open chromatin regions") + # (Only TF motifs with FDR<0.05 are shown)
  geom_hline(aes(yintercept=-log10(0.05), linetype="FDR = 0.05"), color = "black", size=0.3) +
  ylab("TF enrichment\n[-log10(FDR)]") +
  xlab("Transcription factors") +
  labs(fill = "Differential open chromatin regions") +
  theme_bw() +
  theme(text = element_text(size = TextSize),
        axis.text.y = element_text(size = TextSize),
        axis.text.x = element_text(size = TextSize, angle = 45, hjust = 1),
        title = element_text(size = TextSize),
        plot.title = element_text(hjust = 0.5, size = TextSize),
        legend.text = element_text(size = TextSize), 
        legend.title = element_text(size = TextSize),
        legend.margin = margin(0,0,-0.5,0,unit = "lines"),
        legend.key.size = unit(LegendSize, 'lines'),
        legend.position = "top",
        legend.box = "vertical"
  ) +
  scale_linetype_manual(name = "Threshold", values = c("FDR = 0.05" = 2)) + 
  scale_fill_manual(values = c("down-reg in ASD" = "deepskyblue", "up-reg in ASD" = "#F8766D"))

pdf("../../plots/02_TF/2a_01_HOMER_Fig2b.pdf", width = 3.5, height = 2) 
print(Fig2b)
dev.off()

