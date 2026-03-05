rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")
# original directory: results/10_3_8_LinkATACtoGene/v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp_corrected20230227

########## pipeline ########## 
### 1. Identify promoter ATAC (TSS -2kb to +100bp)
### 2. Get all ATAC within consensus promoter Hi-C loops (union of Bulk, NeuNp, NeuNn)
### 3. Get all distal ATAC peaks within 30kb of promoter ATAC peaks
### 4. Calculate ATAC-ATAC cor using RPKM.cqn 
### 5. Build null distribution and calculate p-values
### 6. Save results & generate TableS2
##############################

### 1. Identify promoter ATAC (TSS -2kb to +100bp)
lnames = load("1a_01_CQN.rda") # "RPKM.cqn"   "cqn.fit"    "Covariates" "peak_info"  "counts_tbl"
rm(list = setdiff(ls(), "RPKM.cqn"))

## Load brain expressed genes from Gandal 2022 Nature
library(readxl)
Jill_DEG = readxl::read_excel("../../data/DEG_Gandal2022Nature.xlsx") # original file: ~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx

## Find gene TSS
# BiocManager::install("biomaRt")
library(biomaRt)
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="grch37.ensembl.org") 
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","ensembl_transcript_id",
             "transcript_start", "transcript_end", "strand","transcript_biotype")
geneAnno <- getBM(attributes = getinfo,filters = c("ensembl_gene_id"), values = Jill_DEG$ensembl_gene_id, mart = mart)
unique(geneAnno$strand) # -1 1
geneAnno$TSS = ifelse(geneAnno$strand == "1", geneAnno$transcript_start, geneAnno$transcript_end)
geneAnno$TSS_bin = floor(geneAnno$TSS/10e3) * 10e3 + 5e3
geneAnno$chromosome_name = paste0("chr", geneAnno$chromosome_name)
colnames(geneAnno)[3] = "chr"
geneAnno$strand = ifelse(geneAnno$strand == "-1", "-", "+")
geneAnno$start = ifelse(geneAnno$strand == "+", geneAnno$TSS - 2000, geneAnno$TSS - 100)
geneAnno$end = ifelse(geneAnno$strand == "+", geneAnno$TSS + 100, geneAnno$TSS + 2000)

library(GenomicRanges)
promoter_gr = makeGRangesFromDataFrame(geneAnno)
save(promoter_gr, geneAnno, file = "1c_01_promoter_gr.rda") # original file name: 10_3_8_v5_promoter_gr.rda

library(stringr)
ATAC_location = as.data.frame(str_split_fixed(rownames(RPKM.cqn), "_", 3))
colnames(ATAC_location) = c("chr", "start", "end")
ATAC_location$start = as.integer(ATAC_location$start); ATAC_location$end = as.integer(ATAC_location$end)
ATAC_location$ATACpeak = rownames(RPKM.cqn)
ATAC_gr = makeGRangesFromDataFrame(ATAC_location)

promATACs = as.data.frame(findOverlaps(promoter_gr, ATAC_gr))
promATACs$promoter_ATAC = ATAC_location$ATACpeak[promATACs$subjectHits]
promATACs$ensembl_gene_id = geneAnno$ensembl_gene_id[promATACs$queryHits]
promATACs$chr = geneAnno$chr[promATACs$queryHits]
promATACs$TSS_bin = geneAnno$TSS_bin[promATACs$queryHits]
promATACs = unique(promATACs[,-c(1:2)]) # 29527 promoter ATAC-GE pairs

save(promATACs, file = "1c_01_promATAC.rda") # original file name: promATAC_TSSminus100bptoplus2kb.rda

### 2. Get all ATAC within consensus promoter Hi-C loops (union of Bulk, NeuNp, NeuNn)
# Run 05_xxx.R to get Hi-C loops
lnames = load("../05_HiC/Consensus_promoter_loops_logCPM.rda") # "Bulk_logCPM"  "NeuNp_logCPM" "NeuNn_logCPM" # Original file: ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/2nd_batch_20samples/data_analysis/results/25_nlme_DiffLoopAnalysis/Consensus_promoter_loops_8Bulk7NeuNn9NeuNp_LoopBySample_logCPM_includingOutlier.rda

union_loops = unique(c(rownames(Bulk_logCPM), rownames(NeuNp_logCPM), rownames(NeuNn_logCPM)))
union_loops_df1 = as.data.frame(str_split_fixed(union_loops, "_", 3)) # 67092 loops
union_loops_df1$V2 = as.integer(union_loops_df1$V2); union_loops_df1$V3 = as.integer(union_loops_df1$V3)
union_loops_df1$Loop = union_loops

union_loops_df2 = union_loops_df1
colnames(union_loops_df1)[1:3] = c("chr", "TSS_bin", "Distal_bin")
colnames(union_loops_df2)[1:3] = c("chr", "Distal_bin", "TSS_bin")

promATAC_Loop1 = left_join(promATACs, union_loops_df1)
promATAC_Loop2 = left_join(promATACs, union_loops_df2)
promATAC_Loop1 = promATAC_Loop1[complete.cases(promATAC_Loop1),]
promATAC_Loop2 = promATAC_Loop2[complete.cases(promATAC_Loop2),]
promATAC_Loop = unique(rbind(promATAC_Loop1, promATAC_Loop2)) # 111856 rows

Distal_bin = unique(data_frame(chr = promATAC_Loop$chr, start = promATAC_Loop$Distal_bin - 5e3, end = promATAC_Loop$Distal_bin + 5e3)) # 32511
Distal_bin_gr = annoDF2GR(Distal_bin)
Distal_bin$distal_bin = paste0(Distal_bin$chr, "_", Distal_bin$start + 5e3)
promATAC_Loop$distal_bin = paste0(promATAC_Loop$chr, "_", promATAC_Loop$Distal_bin)

hits = as.data.frame(findOverlaps(Distal_bin_gr, ATAC_gr))
hits$distal_bin = Distal_bin$distal_bin[hits$queryHits]
hits$distal_ATAC = ATAC_location$ATACpeak[hits$subjectHits]
hits = unique(hits[,3:4])

promATAC_Loop = left_join(promATAC_Loop, hits)
promATAC_Loop = promATAC_Loop[complete.cases(promATAC_Loop),] # 154403 promoter-distal ATAC pairs in promoter Hi-C loops
all(promATAC_Loop$promoter_ATAC %in% promATACs$promoter_ATAC) # T

Cor_Gene_ATAC = promATAC_Loop
save(Cor_Gene_ATAC, file = "1c_01_ATACpairsInUnionPL.rda") # original file name: ATACpairsInUnionPL.rda

### 3. Get all distal ATAC peaks within 30kb of promoter ATAC peaks
Cor_Gene_ATAC_within30kb = matrix(nrow = 0, ncol = 3) # ensembl_gene_id, promoter_ATAC, distal_ATAC, cor (omit tss column to reduce calc burden)

for (g in Jill_DEG$ensembl_gene_id) {
  print(g)
  chr = paste0("chr", Jill_DEG$chromosome_name[Jill_DEG$ensembl_gene_id == g])
  
  gene_info = geneAnno[geneAnno$ensembl_gene_id == g,]
  
  # TSS (-100bp to + 2kb)
  tss = gene_info$TSS

  # promoter ATAC
  promATACs_g = unique(promATACs$promoter_ATAC[promATACs$ensembl_gene_id == g])
  
  # TSS +- 30kb
  tss_gr = makeGRangesFromDataFrame(data_frame(chr = chr, start = tss - 30e3, end = tss + 30e3))
  
  # All possible distal ATAC
  hits = as.data.frame(findOverlaps(tss_gr, ATAC_gr))
  distATACs = ATAC_location$ATACpeak[unique(hits$subjectHits)]
  
  # add to data frame
  cur_df = data_frame(ensembl_gene_id = g, promoter_ATAC = rep(promATACs_g, length(distATACs)), distal_ATAC = rep(distATACs, each = length(promATACs_g))) 
  idx_rm = which(cur_df$promoter_ATAC == cur_df$distal_ATAC)
  cur_df = cur_df[-idx_rm, ]
  Cor_Gene_ATAC_within30kb = rbind(Cor_Gene_ATAC_within30kb, cur_df)
}
# Run for 15-30 min
dim(Cor_Gene_ATAC_within30kb) # 221205 3
save(Cor_Gene_ATAC_within30kb, file = "1c_01_ATACpairs_within30kb.rda") # original file: ../v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp_corrected20230227/ATACpairs_within30kb.rda
all(Cor_Gene_ATAC_within30kb$promoter_ATAC %in% promATACs$promoter_ATAC) # T, previously F before the correction

# Combine ATAC pairs in PL and within 30kb
tmp = unique(Cor_Gene_ATAC[,c("ensembl_gene_id", "promoter_ATAC", "distal_ATAC")])
Cor_Gene_ATAC_within30kb = anti_join(Cor_Gene_ATAC_within30kb, tmp) # 190302 3

Cor_Gene_ATAC_within30kb$chr = unname(sapply(Cor_Gene_ATAC_within30kb$promoter_ATAC, function(x) unlist(str_split(x, "_"))[1]))
Cor_Gene_ATAC_within30kb$distal_bin = Cor_Gene_ATAC_within30kb$Distal_bin = Cor_Gene_ATAC_within30kb$TSS_bin = 0
Cor_Gene_ATAC_within30kb$Loop = "Within30kb" # the ones with Hi-C loop is labeled as Loop
Cor_Gene_ATAC_within30kb = Cor_Gene_ATAC_within30kb[,match(colnames(Cor_Gene_ATAC), colnames(Cor_Gene_ATAC_within30kb))]

Cor_Gene_ATAC = rbind(Cor_Gene_ATAC, Cor_Gene_ATAC_within30kb)
save(Cor_Gene_ATAC, file = "1c_01_ATACpairs_PLandWithin30kb.rda") # original file name: ATACpairs_PLandWithin30kb.rda
dim(Cor_Gene_ATAC) # 344705 8

### 4. Calculate ATAC-ATAC cor using RPKM.cqn 
Cor_Gene_ATAC$cor = 0
for (i in 1:nrow(Cor_Gene_ATAC)) {
  if (i %% 20e3 == 0) {print(i)}
  p = Cor_Gene_ATAC$promoter_ATAC[i]
  d = Cor_Gene_ATAC$distal_ATAC[i]
  p_atac = RPKM.cqn[rownames(RPKM.cqn) == p,]
  d_atac = RPKM.cqn[rownames(RPKM.cqn) == d,]
  Cor_Gene_ATAC$cor[i] = cor(p_atac,d_atac)
}
save(Cor_Gene_ATAC, file = "1c_01_ATACpairs_cor.rda") # original file name: Cor_Gene_ATAC_calcedCor.rda
# Run for 15-20 min

# Take a look at the distribution of the correlation between ATAC pairs
hist(Cor_Gene_ATAC$cor) # Observation: It's a skewed distribution, skewed to the right.
mean(Cor_Gene_ATAC$cor, na.rm = T) # PL and within 30kb: 0.23, higher than all promoter-distal ATAC pairs within 500Mb (0.15). 
range(Cor_Gene_ATAC$cor, na.rm = T) # -0.89 to 0.99

### 5. Build null distribution and calculate p-values
# using 10e3 trans ATAC-ATAC cor for each promoter ATAC
uniq_promATAC = unique(Cor_Gene_ATAC$promoter_ATAC); length(uniq_promATAC) # 23768 promoter ATAC peaks

tmp = as.data.frame(table(ATAC_location$chr))
nrow(ATAC_location) - max(tmp$Freq) # 117001, at max can randomly draw 117k trans ATAC peaks.

Cor_Gene_ATAC$p_perPromATAC = Cor_Gene_ATAC$trans_mean_perPromATAC = NA

set.seed(137)
for (i in 1:length(uniq_promATAC)) {
  if (i %% 5e3 == 0) {print(i)}
  p = uniq_promATAC[i]
  chr = unlist(str_split(p, "_"))[1]
  ATAC_chr_pool = which(ATAC_location$chr != chr)
  idx_sample = sample(ATAC_chr_pool, 10e3) # 20e3 doesn't make much difference
  null_atac = ATAC_location$ATACpeak[idx_sample]
  
  null_RPKM = t(RPKM.cqn[rownames(RPKM.cqn) %in% null_atac,])
  p_RPKM = RPKM.cqn[rownames(RPKM.cqn) == p,]
  null_cor = cor(p_RPKM, null_RPKM)
  # hist(null_cor) # It's also skewed, and sometimes irregular shaped.
  
  m = mean(null_cor) # 0.04 for chrX_99890783_99892116

  idx = which(Cor_Gene_ATAC$promoter_ATAC == p)
  Cor_Gene_ATAC$trans_mean_perPromATAC[idx] = m
  Cor_Gene_ATAC$p_perPromATAC[idx] = sapply(Cor_Gene_ATAC$cor[idx], function(x) length(which(null_cor > x))/10e3) # one-tail test
}
# Run for 5 min

hist(Cor_Gene_ATAC$p_perPromATAC, breaks = 50) # wow, very enriched for p< 0.05
Cor_Gene_ATAC$fdr_perPromATAC = p.adjust(Cor_Gene_ATAC$p_perPromATAC, method = "fdr")
save(Cor_Gene_ATAC, file = "1c_01_ATACpairs_corP.rda") # original file name: Cor_Gene_ATAC_calcedCorP.rda

# Use FDR < 0.1 for significant promoter-distal ATAC pairs
sig_CorATAC = Cor_Gene_ATAC[Cor_Gene_ATAC$fdr_perPromATAC < 0.1,] # 25184 ATAC-ATAC pairs

# Attach differential ATAC results
lnames = load("1a_03_DiffATAC.rda")
sig_CorATAC$promoterATAC_logFC = DiffATAC$ATAC_logFC[match(sig_CorATAC$promoter_ATAC, rownames(DiffATAC))]
sig_CorATAC$promoterATAC_FDR = DiffATAC$ATAC_FDR[match(sig_CorATAC$promoter_ATAC, rownames(DiffATAC))]
sig_CorATAC$distalATAC_logFC = DiffATAC$ATAC_logFC[match(sig_CorATAC$distal_ATAC, rownames(DiffATAC))]
sig_CorATAC$distalATAC_FDR = DiffATAC$ATAC_FDR[match(sig_CorATAC$distal_ATAC, rownames(DiffATAC))]
sig_CorATAC$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(sig_CorATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
sig_CorATAC$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(sig_CorATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

sig_CorATAC_uniq = unique(sig_CorATAC[,c("distal_ATAC", "ensembl_gene_id", "distalATAC_logFC", "distalATAC_FDR", "WholeCortex_ASD_logFC", "WholeCortex_ASD_FDR")]) # 21067 rows

### 6. Save results & generate TableS2
# Attach differential ATAC results to promoter ATAC peaks
promATAC = promATACs[,1:2]
colnames(promATAC)[1] = "ATACpeak"
promATAC$ATAC_logFC = DiffATAC$ATAC_logFC[match(promATAC$ATACpeak, rownames(DiffATAC))]
promATAC$ATAC_FDR = DiffATAC$ATAC_FDR[match(promATAC$ATACpeak, rownames(DiffATAC))]
promATAC$WholeCortex_ASD_logFC = Jill_DEG$WholeCortex_ASD_logFC[match(promATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]
promATAC$WholeCortex_ASD_FDR = Jill_DEG$WholeCortex_ASD_FDR[match(promATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

save(sig_CorATAC, sig_CorATAC_uniq, promATAC, file = "1c_01_CHART_assigned_distalATAC_gene_pairs.rda") # original file name: v5_ATACatacCor_PL_PromoterAsTSSminus2kbplus100bp/Distal_And_Promoter_ATAC_GE_pairs_FDRcorATAC01.rda

## Generate Supplementary Table 2
rm(list = ls())
lnames = load("1c_01_CHART_assigned_distalATAC_gene_pairs.rda")
length(unique(sig_CorATAC$ensembl_gene_id)) # 8865
length(unique(sig_CorATAC$distal_ATAC)) # 15420

sig_CorATAC = sig_CorATAC[,c("ensembl_gene_id", "promoter_ATAC", "distal_ATAC", "Loop", "cor", "fdr_perPromATAC")] # Loop contains either the promoter loop location info or "Within30kb"
colnames(sig_CorATAC)[c(4,6)] = c("loop", "FDR")

# Add gene_name
Jill_DEG = readxl::read_excel("../../data/DEG_Gandal2022Nature.xlsx") # original file: ~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx
sig_CorATAC$gene_name = Jill_DEG$external_gene_name[match(sig_CorATAC$ensembl_gene_id, Jill_DEG$ensembl_gene_id)]

sig_CorATAC = sig_CorATAC[,c(7, 1:6)]

# sort the table by promoter ATAC start and end positions
sig_CorATAC$chr = 
  sapply(sig_CorATAC$promoter_ATAC, function(x) substr(unlist(strsplit(x, "_", 3))[1], 4, 6))
sig_CorATAC$promoter_ATAC_start = 
  sapply(sig_CorATAC$promoter_ATAC, function(x) as.integer(unlist(strsplit(x, "_", 3))[2]))
sig_CorATAC$promoter_ATAC_end = 
  sapply(sig_CorATAC$promoter_ATAC, function(x) as.integer(unlist(strsplit(x, "_", 3))[3]))
sig_CorATAC$chr[sig_CorATAC$chr == "X"] = 23
sig_CorATAC$chr[sig_CorATAC$chr == "Y"] = 24
sig_CorATAC$chr = as.integer(sig_CorATAC$chr)

sig_CorATAC = sig_CorATAC %>%
  arrange(chr, promoter_ATAC_start, promoter_ATAC_end)
unique(sig_CorATAC$chr) # looks good

TableS2 = sig_CorATAC[,1:7]

# Write to Excel
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, '1c_01_CHART_TableS2')
header_st <- createStyle(textDecoration = "Bold")
writeData(wb, '1c_01_CHART_TableS2', TableS2, headerStyle = header_st)
saveWorkbook(wb, '1c_01_CHART_TableS2.xlsx', overwrite = TRUE)
# Original file: /Users/dhglab/Desktop/ASD_manuscript/02_Submission/06_SupplementaryTables_v14/Supplementary Table 2.xlsx


