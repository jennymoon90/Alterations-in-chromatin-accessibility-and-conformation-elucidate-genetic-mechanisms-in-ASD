rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")

########## pipeline ########## 
### 1. Load CHART-assigned DAR target genes
### 2. GO analysis on up/down-reg DARs at promoter and distal regions
### 3. Bubble plot
##############################

### 1. Load CHART-assigned DAR target genes
lnames = load("1c_01_CHART_assigned_distalATAC_gene_pairs.rda") # "sig_CorATAC"      "sig_CorATAC_uniq" "promATAC"
# original file: ../v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda

### 2. GO analysis on up/down-reg DARs at promoter and distal regions

## Get the target gene sets for each category of interest
any(is.na(promATAC$ensembl_gene_id)) # F
any(is.na(sig_CorATAC_uniq$ensembl_gene_id)) # F

Promoter_down_geneid = unique(promATAC$ensembl_gene_id[which(promATAC$ATAC_logFC < 0 & promATAC$ATAC_FDR < 0.05)]) # 533 genes
Promoter_up_geneid = unique(promATAC$ensembl_gene_id[which(promATAC$ATAC_logFC > 0 & promATAC$ATAC_FDR < 0.05)]) # 960 genes
DistalA_down_geneid = unique(sig_CorATAC_uniq$ensembl_gene_id[which(sig_CorATAC_uniq$distalATAC_logFC < 0 & sig_CorATAC_uniq$distalATAC_FDR < 0.05)]) # 791 genes
DistalA_up_geneid = unique(sig_CorATAC_uniq$ensembl_gene_id[which(sig_CorATAC_uniq$distalATAC_logFC > 0 & sig_CorATAC_uniq$distalATAC_FDR < 0.05)]) # 511 genes

## Use all brain-expressed genes as background
library(readxl)
Jill_DEG = readxl::read_excel("../../data/DEG_Gandal2022Nature.xlsx") # original file: ~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx

## GO analysis and plot using Level 4 terms
library(clusterProfiler) # T Wu etc. The Innovation. 2021
gost_PromoterDown_genes = enrichGO(gene = Promoter_down_geneid, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", ont = "BP", universe = Jill_DEG$ensembl_gene_id)
gost_PromoterDown_genes_filtered = gofilter(gost_PromoterDown_genes, level = c(4,5))
GO_PromoterDown_genes = gost_PromoterDown_genes_filtered@result
GO_PromoterDown_genes$minusLog10Pval = -log10(GO_PromoterDown_genes$pvalue) 
GO_PromoterDown_genes = GO_PromoterDown_genes[GO_PromoterDown_genes$qvalue < 0.15,] # 0 term

gost_PromoterUp_genes = enrichGO(gene = Promoter_up_geneid, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", ont = "BP", universe = Jill_DEG$ensembl_gene_id)
GO_PromoterUp_genes = gost_PromoterUp_genes@result
gost_PromoterUp_genes_filtered = gofilter(gost_PromoterUp_genes, level = c(4,5))
GO_PromoterUp_genes = gost_PromoterUp_genes_filtered@result
GO_PromoterUp_genes$minusLog10Pval = -log10(GO_PromoterUp_genes$pvalue) 
GO_PromoterUp_genes = GO_PromoterUp_genes[GO_PromoterUp_genes$qvalue < 0.15,] # 1 term

gost_DistalDown_genes = enrichGO(gene = DistalA_down_geneid, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", ont = "BP", universe = Jill_DEG$ensembl_gene_id)
gost_DistalDown_genes_filtered = gofilter(gost_DistalDown_genes, level = c(4,5))
GO_DistalDown_genes = gost_DistalDown_genes_filtered@result
GO_DistalDown_genes$minusLog10Pval = -log10(GO_DistalDown_genes$pvalue) 
GO_DistalDown_genes = GO_DistalDown_genes[GO_DistalDown_genes$p.adjust < 0.05,] # 4 terms

gost_DistalUp_genes = enrichGO(gene = DistalA_up_geneid, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", ont = "BP", universe = Jill_DEG$ensembl_gene_id)
gost_DistalUp_genes_filtered = gofilter(gost_DistalUp_genes, level = c(4,5))
GO_DistalUp_genes = gost_DistalUp_genes_filtered@result
GO_DistalUp_genes$minusLog10Pval = -log10(GO_DistalUp_genes$pvalue) 
GO_DistalUp_genes = GO_DistalUp_genes[GO_DistalUp_genes$qvalue < 0.1,] # 5 terms

## Organize the dataframe for plot
GO_PromoterDown_genes$ATAC_type = "Promoter_Down"
GO_PromoterUp_genes$ATAC_type = "Promoter_Up"
GO_DistalDown_genes$ATAC_type = "Distal_Down"
GO_DistalUp_genes$ATAC_type = "Distal_Up"

GO_df = rbind(GO_PromoterDown_genes, GO_PromoterUp_genes, GO_DistalDown_genes, GO_DistalUp_genes)
length(unique(GO_df$Description)) # 15 (final decision) unique GO term names

GO_df$minFDR = sapply(GO_df$Description, function(x) min(GO_df$p.adjust[GO_df$Description == x]))
GO_df$ATAC_type = factor(GO_df$ATAC_type, levels = c("Promoter_Up", "Promoter_Down", "Distal_Up", "Distal_Down"))
GO_df_selectedTerms = GO_df %>%
  arrange(ATAC_type, desc(minusLog10Pval))

## Shorten the GO term name
idx = which(nchar(GO_df_selectedTerms$Description) > 48) # 8
lengths(strsplit(GO_df_selectedTerms$Description[idx], ' ')) # 6 words
library(tidyverse)
word_split <- function(x, side="left", sep=" ") {
  words <- strsplit(as.character(x), sep)
  nwords <- lengths(words)
  if(side=="left") {
    start <- 1
    end <- ceiling(nwords/2)
  } else if (side=="right") {
    start <- ceiling((nwords+1)/2)
    end <- nwords
  }
  cw <- function(words, start, stop) paste(words[start:stop], collapse=sep)
  pmap_chr(list(words, start, end), cw)
}
left_words = word_split(GO_df_selectedTerms$Description[idx])
right_words = word_split(GO_df_selectedTerms$Description[idx], side = "right")
GO_df_selectedTerms$Description[idx] = paste0(left_words, "\n", right_words)

GO_df_selectedTerms$Description = factor(GO_df_selectedTerms$Description, levels = rev(unique(GO_df_selectedTerms$Description)))

save(GO_df, GO_df_selectedTerms, file = "1c_03_DAR_GOterms.rda") # original file: $DIR/10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp_DistalcorWithPromoterTSSplusminus2kb/DiffATAC_PromoterOrDistalAssigned_geneGOenrichment_L45terms.rda

### 3. Bubble plot
rm(list = ls())
lnames = load("1c_03_DAR_GOterms.rda")

GO_df_selectedTerms$minusLog10FDR = -log10(GO_df_selectedTerms$qvalue)
# In order to have the "Promoter_Down" category show up in legend with a color, add the NA line of this category
GO_df_selectedTerms[nrow(GO_df_selectedTerms) + 1,] = NA
GO_df_selectedTerms$ATAC_type[nrow(GO_df_selectedTerms)] = "Promoter_Down"

TextSize = 5
LegendSize = 0.5

Fig1h_bubble = 
  GO_df_selectedTerms[! is.na(GO_df_selectedTerms$Description),] %>%
  ggplot() +
  geom_point(aes(x = Description, y = ATAC_type, col = ATAC_type, size = minusLog10FDR)) + 
  scale_size(range = c(0.5,3.5)) + # adjust point size to fit the figure size! 
  theme_bw() +
  coord_flip() + 
  scale_linetype_manual(name = "Theshold", values = 2) +
  ylab("DAR type") +
  labs(size = "-log10(FDR)", col = "DAR type") +
  scale_color_manual(labels = c("Promoter Up", "Promoter Down", "Distal Up", "Distal Down"), values = c("#F8766D", "darkgoldenrod1", "#00BFC4", "#C77CFF"), drop = F) + #\n(linked to promoter)
  ggtitle("GO enrichment of target genes of DARs") +
  theme(plot.title = element_text(hjust = 0.65, size = TextSize),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = TextSize),
        text = element_text(size = TextSize),
        title = element_text(size = TextSize),
        legend.text = element_text(size = TextSize),
        legend.key.size = unit(LegendSize, 'lines'),
        legend.margin = margin(0,unit = "lines"),
        legend.spacing.y = unit(LegendSize/2, 'lines'),
        panel.grid.minor = element_line(linewidth = 0.25), 
        panel.grid.major = element_line(linewidth = 0.5) #,
  )

pdf("../../plots/01_DAR/1c_03_GO_DARs_Fig1h.pdf", width = 3, height = 2) # original file: DiffATAC_PromoterOrDistalAssigned_geneGOenrichment_L45terms.pdf
Fig1h_bubble
dev.off()


