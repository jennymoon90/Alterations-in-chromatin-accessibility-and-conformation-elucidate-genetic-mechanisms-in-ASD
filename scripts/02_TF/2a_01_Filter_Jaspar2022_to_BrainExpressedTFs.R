rm(list=ls())
options(stringsAsFactors=F)

library(readxl)
library(stringr)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/utils/")
# original directory: ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/1st_batch_parietal/4_Run_for_all_samples_goal_10kb/results/

########## pipeline ########## 
### 1. Download JAPSAR2022 vertebrate motifs
### 2. Only keep TFs that are expressed in human brain
### 3. Call consensus motif sequence from frequency matrix
### 4. Calculate Log odds detection threshold for each motif
### 5. Save the resulting .motifs file
##############################

### 1. Download JAPSAR2022 vertebrate motifs 
# from https://jaspar.genereg.net/downloads/ -> Vertebrates tab -> single batch file (txt) -> PFMs(non-redundant) JASPAR -> save link as utils/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt # original file: ~/Documents/Documents/Geschwind_lab/LAB/Computational_tools/HOMER/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt
JASPAR = read.table("JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt", col.names = paste0("V",seq_len(40)), fill = T) # original file: ~/Documents/Documents/Geschwind_lab/LAB/Computational_tools/HOMER/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt
names_idx = grep(pattern = ">", x = JASPAR$V1)
(tmp = (names_idx -1)/5) # check that they are all integers, then pick the smallest N for seq_len(N) in JASPAR = read.table()

### 2. Only keep TFs that are expressed in human brain
DEG = read_excel("../data/DEG_Gandal2022Nature.xlsx") # 24836 genes in DEG data from Jill's pan-cortical study. # original file: ~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/logFC_SupplementaryTable3.xlsx
Br_expressed_genes = DEG$external_gene_name

JASPAR_TF_split = as.data.frame(str_split_fixed(JASPAR$V2, "::", 2))
idx1 = which(JASPAR_TF_split$V1 %in% Br_expressed_genes)
idx2 = which(JASPAR_TF_split$V2 %in% Br_expressed_genes)
idx = unique(c(idx1,idx2)) # 556 brain expressed transcription factors in JASPAR2022 database.

idx_withMatrix = c(idx, idx+1, idx+2,idx+3,idx+4) # the 4 rows after idx are frequency matrix, we need them
idx_withMatrix = sort(idx_withMatrix)
JASPAR = JASPAR[idx_withMatrix,]

### 3. Call consensus motif sequence from frequency matrix
idx = grep(">",JASPAR$V1)
for (i in idx) {
  pfm = JASPAR[(i+1):(i+4),]
  end = which(pfm[1,] == "]")
  consensus = rep("N", end-3)
  
  for (j in 3:(end-1)) {
    # MEME-consensus algorithm, refer to http://meme-suite.org/doc/motif-consensus.html
    column = pfm[,c(1,j)]
    column[,2] = as.numeric(column[,2])
    column$freq = column[,2]/sum(column[,2])
    column = column[order(column$freq, decreasing = T),]
    column = column[column[,2] >= column[1,2]/2,] # column[1,2] represents the maximum frequency
    consensus[j-2] = 
      ifelse(nrow(column) == 1, column[1,1],
             ifelse(sum(column$freq) >= 0.5, "N", column[1,1]))
  }
  
  consensus # looks concordant to JASPAR motif sequence logo.
  JASPAR$V3[i] = paste(consensus, collapse = '')
}

### 4. Calculate Log odds detection threshold for each motif
TFnames = JASPAR[idx,1:3]
colnames(TFnames) = c("JASPAR_motif_ID","TF_name","consensus_sequence")
TFnames$JASPAR_motif_ID = substr(TFnames$JASPAR_motif_ID,2,nchar(TFnames$JASPAR_motif_ID))
TFnames$match_jaspar.motifs = paste0(TFnames$TF_name,"/",TFnames$JASPAR_motif_ID,"/Jaspar")
TFnames$code = paste("seq2profile.pl", TFnames$consensus_sequence, 0, TFnames$match_jaspar.motifs, ">",TFnames$JASPAR_motif_ID)
TFnames$code = paste0(TFnames$code,".motif")
TFnames$code = gsub("\\(var.2\\)","-var2",TFnames$code)

write.table(TFnames$code, file = "../scripts/02_TF/2a_02_HOMER_seq2profile_Jaspar2022_ind.sh", quote = F, row.names = F, col.names = F, sep = "\t") # original file ~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_ASD/Hi-C/1st_batch_parietal/4_Run_for_all_samples_goal_10kb/results/14_GETPM_Nloops/motifs/HOMER_seq2profile_Jaspar2022.sh

## Go to terminal -> follow the steps in 2a_03_HOMER_seq2profile_Jaspar2022_all.sh to run 2a_02_HOMER_seq2profile_Jaspar2022_ind.sh and cat all .motif files to one file
# Originally: Go to terminal -> follow the steps in 14_11_HOMER_scan_motif_in_promoterATAC_4categories_Jaspar2022.sh to run HOMER_seq2profile_Jaspar2022.sh and cat all .motif files to one file

LogOdds = read.table("res_seq2profile_Jaspar2022/Jaspar2022_BrExp_withHOMERlogodds.motif", header = F, sep = "\t", fill = T)
LogOdds$V2 = gsub("-var2","\\(var.2\\)",LogOdds$V2)

## Correct the frequency matrix with JASPAR outputs
idx = grep(">",LogOdds$V1)
for (i in idx) {
  motif_id = str_split_fixed(LogOdds$V2[i],"/",3)[,2]
  JASPAR_rownumber = grep(motif_id,JASPAR$V1)
  pfm = JASPAR[(JASPAR_rownumber+1):(JASPAR_rownumber+4),]
  end = which(pfm[1,] == "]")
  pfm = pfm[,3:(end-1)] # rows are A,C,G,T in order
  pfm = t(pfm)
  pfm = as.data.frame(apply(pfm,2,as.numeric))
  st = apply(pfm,1,sum)
  pfm = sweep(pfm,1,st,FUN="/")
  pfm = round(pfm,digits = 3)
  
  LogOdds[(i+1):(i+nrow(pfm)),1:4] = pfm
}

LogOdds$V4[is.na(LogOdds$V4)] = ""

### 5. Save the resulting .motifs file
write.table(LogOdds,file = "jaspar2022_BrExpTF.motifs", quote = F, row.names = F, col.names = F, sep = "\t")

