### Goal: Get top PCs that explain total 80% variance in the CQN normalized ATAC peak count, identify biological and technical covariates that significantly correlate with these top PCs.

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)

setwd("~/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/")

########## pipeline ########## 
### 1. Load data
  ## 1) CQN normalized peak count data
  ## 2) ATAC library biological and technical covariates
### 2. Remove technical replicates for building linear model
### 3. Covariate selection
  ## 1) Identify top PCs
  ## 2) Correlation of covariates with the top PCs
  ## 3) Find out which Covariates significantly correlate with the topPCs
  ## 4) Select covariates
### 4. Run linear regression (use linear mixed model, including random effect of subjects)
### 5. Save differential ATAC results (ASD vs. CTL)
##############################

### 1. Load data
## 1) CQN normalized peak count data
lnames = load("1a_01_CQN.rda") # original file name: CQN.rda
rm(list = setdiff(ls(), "RPKM.cqn"))

## 2) ATAC library biological and technical covariates
lnames = load("../../data/ATAC_library_covariates.rda") # df_wAllCov. original file name: 3_QC/BiolTechCov_and_QC.rda
Covariates = df_wAllCov; rm(df_wAllCov)
Covariates = Covariates[match(colnames(RPKM.cqn), Covariates$sample_name),] # 38 rows/samples

idx_rm = which(grepl("_quality", colnames(Covariates)))
Covariates = Covariates[,-idx_rm]
idx_rm = which(colnames(Covariates) %in% c("QC_summary", "Subject_ID", "Subject_ID_2", "sample_name"))
Covariates = Covariates[,-idx_rm]

Covariates$PMI = as.numeric(Covariates$PMI) 
Covariates$PMI[is.na(Covariates$PMI)] = mean(Covariates$PMI, na.rm = T)
Covariates$Diagnosis = factor(Covariates$Diagnosis, levels = c("CTL", "ASD"))
Covariates$Sex = factor(Covariates$Sex, levels = c("M", "F"))
Covariates$Cortex = factor(Covariates$Cortex, levels = c("Parietal", "Frontal", "Temporal"))
Covariates$Age2 = Covariates$Age^2

colnames(Covariates) = gsub("_", ".", colnames(Covariates))

# Limit biol and tech cov to the relevant ones selected based on the "unorganized" scripts
colnames(Covariates)[which(colnames(Covariates) == "subject")] = "Subject"
Covariates = Covariates[,c("Subject","Diagnosis", "Age", "Age2", "PMI", "Sex", "BrainBank", "Batch", "Cortex", "Region", "tssenrich.score", "FRiP", "picard.alignment.PCT.READS.ALIGNED.IN.PAIRS", "picard.duplication.PERCENT.DUPLICATION")]
colnames(Covariates)[(ncol(Covariates)-1):ncol(Covariates)] = c("PCT.READS.ALIGNED.IN.PAIRS", "PERCENT.DUPLICATION")

### 2. Remove technical replicates for building linear model
# But when I run lme, I shall include Subject as a random effect
TechRep_rm = c("B4334-1", "B4721A", "B5000B", "B5718D", "B5813B", "CQ56-2")
RPKM.cqn2 = RPKM.cqn[,which(! colnames(RPKM.cqn) %in% TechRep_rm)]
Covariates2 = Covariates[which(! rownames(Covariates) %in% TechRep_rm),]

stopifnot(rownames(Covariates) == colnames(RPKM.cqn))
stopifnot(rownames(Covariates2) == colnames(RPKM.cqn2))

### 3. Covariate selection
## 1) Identify top PCs
norm <- t(scale(t(RPKM.cqn2),scale=F))
PC <- prcomp(norm,center=FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[c(1:14)]) ## Top 14 PCs explain 79.7% of total variance.
topPC <- PC$rotation[,1:14]

## 2) Correlation of covariates with the top PCs
mod_mat_expr = paste(c(colnames(Covariates2)[-1]), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates2)[,-1]
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
idx_spearman = ncol(topPC) + c(1, 5:15)
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

## 3) Find out which Covariates significantly correlate with the topPCs
Cor_sig = matrix(nrow = nrow(Cor), ncol = ncol(Cor))
rownames(Cor_sig) = colnames(Cor_sig) = colnames(Cor)
for (i in 1:ncol(mod_mat_withPC)) {
  for (j in 1:ncol(mod_mat_withPC)) {
    tmp = cor.test(mod_mat_withPC[,i], mod_mat_withPC[,j])
    Cor_sig[i,j] = tmp$p.value
    # if (tmp$p.value < 0.05) {print(paste(i, colnames(Cor)[i], j, colnames(Cor)[j]))}
  }
}

n_tests = (ncol(Cor_sig)-1)^2 - (ncol(topPC)-1)^2
(p_cor_threshold = 0.05/n_tests)
Cov_pool_idx = which(apply(Cor_sig[1:ncol(topPC),], 2, function(x) any(x < p_cor_threshold)))
(Cov_pool = colnames(Cor_sig)[Cov_pool_idx]) # FRiP
# Unfortunately diagnosis does not significantly correlate with any topPCs, but this is fine.

# Maybe bonforroni correction of the p-val is too stringent, try use fdr and see how many covaraites are candidates
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}

Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # No DiagnosisASD! Age, Age2, Batch, CortexFrontal, CortexTemporal, RegionBA38, RegionBA44_45, RegionBA9, tssenrich.score, FRiP, picard.alignment.PCT.READS.ALIGNED.IN.PAIRS, picard.duplication.PERCENT.DUPLICATION.
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.1)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Still no DiagnosisASD! Age, Age2, BrainBankNICHD-BTB, Batch, CortexFrontal, CortexTemporal = RegionBA38, RegionBA44_45, RegionBA9, tssenrich.score, FRiP, picard.alignment.PCT.READS.ALIGNED.IN.PAIRS, picard.duplication.PERCENT.DUPLICATION.
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.2)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # Now there is DiagnosisASD! Use this fdr<0.2 as threshold. DiagnosisASD, Age, Age2, BrainBankNICHD-BTB, Batch, CortexFrontal, CortexTemporal = RegionBA38, RegionBA44_45, RegionBA9, tssenrich.score, FRiP, picard.alignment.PCT.READS.ALIGNED.IN.PAIRS, picard.duplication.PERCENT.DUPLICATION.

save(Cor, topPC, file = "1a_03_correlation_topPC_covariates.rda") # original file name: Corrplot_topPCofCQNnormedRPKM_Covariates.rda

## Corrplot (the 2nd plot of each denotes fdr<0.2 by *)
library(corrplot)
pdf(paste0("../../plots/01_DAR/1a_03_Corrplot_topPC_covariates_SupFig2.pdf"), height = 35, width = 45) # original file: Corrplot_topPCofCQNnormedRPKM_Covariates.pdf
corrplot(Cor,method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor, p.mat = Cor_fdr, insig = "label_sig", sig.level = 0.2, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2) 
corrplot(Cor[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[1:ncol(topPC), (ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = 0.2, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey",
         tl.cex = 2, cl.cex = 2, number.cex = 2)
corrplot(Cor[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)], p.mat = Cor_fdr[(ncol(topPC) + 1):ncol(Cor),(ncol(topPC) + 1):ncol(Cor)], insig = "label_sig", sig.level = 0.2, pch.col = "white",tl.pos = "lt",tl.col = "black", tl.srt = 45, tl.cex = 2, cl.cex = 2)
dev.off()

## 4) Select covariates
## I know that Batch/Cortex/Region/BrainBank are correlated
unique(Covariates[,c("Cortex", "Region", "Batch")])
# Frontal       BA9         2
# Temporal      BA38        2
# Parietal      BA7         1
# Parietal      BA3_1_2_5   1
# Frontal       BA44_45     2
# Parietal      BA39_40     1
unique(Covariates[,c("BrainBank", "Batch")])
# ABN     2
# NICHD-BTB     2
# NICHD-BTB     1
# Harvard-ATP     1
# Harvard-ATP     2
unique(Covariates[,c("Cortex", "Region", "Batch", "BrainBank")])
# Frontal        BA9     2         ABN
# Temporal      BA38     2   NICHD-BTB
# Parietal       BA7     1   NICHD-BTB
# Parietal BA3_1_2_5     1   NICHD-BTB
# Parietal       BA7     1 Harvard-ATP
# Frontal    BA44_45     2   NICHD-BTB
# Frontal        BA9     2   NICHD-BTB
# Temporal      BA38     2 Harvard-ATP
# Parietal BA3_1_2_5     1 Harvard-ATP
# Parietal   BA39_40     1 Harvard-ATP

Cov_pool # All candidate covariates

## First include all candidate covariates that cor > 0.3 with PC1 or PC2 and that are not significantly correlated with each other.
# Prioritize biological before technical covariate.
test_row = which(rownames(Cor_sig) == "DiagnosisASD")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none.
# Include DiagnosisASD

test_row = which(rownames(Cor_sig) == "Age")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Age2
test_row = which(rownames(Cor_sig) == "Age2")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Age
# Choose between Age and Age2. Age and Age2 are almost 0.98 correlated. Just include Age.

test_row = which(rownames(Cor_sig) == "Batch")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # CortexFrontal, RegionBA9, PCT.READS.ALIGNED.IN.PAIRS, PERCENT.DUPLICATION
test_row = which(rownames(Cor_sig) == "RegionBA9")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # Batch, CortexFrontal, PERCENT.DUPLICATION
test_row = which(rownames(Cor_sig) == "RegionBA38")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # CortexTemporal
# Include RegionBA38. Now I either include Batch or BA9. Batch would take care of the two technical covariates, and the sum of its |cor| with top 2 PCs are higher than BA9, so take in Batch and may spare BA9 if VIF > 2.5.
sum(abs(Cor[1:2,which(colnames(Cor) == "Batch")])) # 0.848
sum(abs(Cor[1:2,which(colnames(Cor) == "RegionBA9")])) # 0.778

test_row = which(rownames(Cor_sig) == "tssenrich.score")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # FRiP
test_row = which(rownames(Cor_sig) == "FRiP")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # PC2, tssenrich.score
# Choose between tssenrich.score and FRiP. Take in tssenrich.score and may spare FRiP if VIF > 2.5
sum(abs(Cor[1:2,which(colnames(Cor) == "tssenrich.score")])) # 1.048
sum(abs(Cor[1:2,which(colnames(Cor) == "FRiP")])) # 1.004

test_row = which(rownames(Cor_sig) == "BrainBankNICHD-BTB")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # None
test_row = which(rownames(Cor_sig) == "RegionBA44_45")
picard_idx = which(Cor_sig[test_row,] < p_cor_threshold)
(picard_cor = colnames(Cor_sig)[picard_idx]) # none

## Check VIF
library(olsrr)
# Base model
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score, data = cur_data)"
# All VIF < 2

# Will BA9 cause a problem?
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score + RegionBA9, data = cur_data)"
# Batch, BA38 and BA9 have VIF > 4. Do not include BA9.

# Will FRiP cause a problem?
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score + FRiP, data = cur_data)"
# All VIF < 2.5. It's fine to take in FRiP

# Will DUPLICATION cause a problem?
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score + PERCENT.DUPLICATION, data = cur_data)"
# Batch and PERCENT.CUPLICATION VIF > 6. Do not inclue DUPLICATION. Same for ALIGNED.
# Now examine the two covariates that have low cor with the top 2 PCs:

# Will RegionBA44_45 or NICHD cause a problem?
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score + FRiP + RegionBA44_45, data = cur_data)"
# All VIF < 2.5. It's fine to take in RegionBA44_45 alone
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score + FRiP + BrainBankNICHD, data = cur_data)"
# All VIF < 2.5. It's fine to take in BrainBankNICHD alone
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + tssenrich.score + FRiP + RegionBA44_45 + BrainBankNICHD, data = cur_data)"
# FRiP VIF > 2.5. Can't take in both RegionBA44_45 and BrainBankNICHD. Decide to take in RegionBA44_45, as it has higher sum(|cor|) with the top 2 PCs.
sum(abs(Cor[1:2,which(colnames(Cor) == "RegionBA44_45")])) # 0.331
sum(abs(Cor[1:2,which(colnames(Cor) == "BrainBankNICHD-BTB")])) # 0.196

# --- run this section for each expression_lm to check VIF ---
i=1
cur_data = RPKM.cqn2[i,]
cur_data = as.data.frame(cbind(cur_data, mod_mat))
colnames(cur_data)[1] = c("logRPKM")
colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # VIF over 5 is warning sign
# ------------------------------------------------------------

## Final model: 
expression_lm = "lm(logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38 + RegionBA44_45 + tssenrich.score + FRiP, data = cur_data)"

### 4. Run linear regression (use linear mixed model, including random effect of subjects)
library(nlme)
runlme <- function(thisdat,expression) {
  lm1 <- eval(parse(text=expression));
  lm1.summary = summary(lm1)
  tabOut <- lm1.summary$coefficients$fixed
  lm1.anova = anova(lm1)
  return(list(tabOut, lm1.anova))
}

Covariates$Subject = factor(Covariates$Subject)
Covariates$RegionBA9 = ifelse(Covariates$Region == "BA9", 1, 0)
Covariates$RegionBA38 = ifelse(Covariates$Region == "BA38", 1, 0)
Covariates$RegionBA44_45 = ifelse(Covariates$Region == "BA44_45", 1, 0)
Covariates$BrainBankNICHD = ifelse(Covariates$BrainBank == "BrainBankNICHD-BTB", 1, 0)

expression_model = "lme(logRPKM ~ Diagnosis + Age + Batch + RegionBA38 + RegionBA44_45 + tssenrich.score + FRiP, random=~1|Subject, data = cur_data)"
n = 7
p = magnitude = matrix(nrow = nrow(RPKM.cqn), ncol = n) 

for (i in 1:nrow(RPKM.cqn)) {
  if (i %% 50000 == 0) {print(paste0("Done ", i, "th loop"))}
  cur_data = RPKM.cqn[i,]
  cur_data = as.data.frame(cbind(cur_data, Covariates))
  colnames(cur_data)[1] = c("logRPKM")
  #colnames(cur_data)[which(colnames(cur_data) == "BrainBankNICHD-BTB")] = "BrainBankNICHD"
  lm1.out <- try(runlme(cur_data,expression_model),silent=F)
  
  if (substr(lm1.out[1],1,5)!="Error") {
    tabOut <- lm1.out[[1]]
    lm1.anova = lm1.out[[2]]
    magnitude[i,] <- tabOut[-1]
    p[i,] <- lm1.anova[-1,"p-value"]
  } else {
    cat('Error in LME of ATAC peak', i, rownames(RPKM.cqn)[i],'\n')
    cat('Setting P-value=NA,Beta value=NA, and SE=NA\n')
    magnitude[i,] <- p[i,] <- NA
  }
}
# Runs for a while (15min) with many Error messages. Fine, cuz this is linear mixed model.
length(which(is.na(p[,1]))) 
# 75 ATAC peaks show error (Error in lme.formula(logRPKM ~ Diagnosis + Age + Batch + RegionBA38 +  : nlminb problem, convergence error code = 1; message = singular convergence (7)). # slightly different from my original analysis (65 error messages), probably due to different package version (re-installed a few).

tabOut
colnames(p) = colnames(magnitude) = names(tabOut)[-1]
rownames(p) = rownames(magnitude) = rownames(RPKM.cqn)
save(p, magnitude, file = "1a_03_lme_p_magnitude_new.rda") # original file name: p_magnitude_lmeDiagnosisASD_ASD_Batch_RegionBA38_RegionBA44.45_tssenrich.score_FRiP.rda

## Check p-value distribution
pdf("../../plots/01_DAR/1a_03_lme_pDiagnosis_histogram.pdf", height = 5, width = 8) # original file name: HistPval_lmeDiagnosisASD_Age_Batch_RegionBA38_RegionBA44.45_tssenrich.score_FRiP.pdf
hist(p[,1], breaks = 50, xlab = "p-value by Diagnosis", main = "Linear mixed model: logRPKM ~ DiagnosisASD + Age + Batch + RegionBA38\n + RegionBA44_45 + tssenrich.score + FRiP + 1|Subject")
dev.off()
# Looks amazing, very sharp p<0.05

## Check if there is any significant hits (fdr < 0.05)
fdr = p.adjust(p[,1], method = "fdr")
range(fdr, na.rm = T)
length(which(fdr < 0.05))
length(which(fdr < 0.1))
length(which(fdr < 0.2))
# FDR range from 0.007 to 1. 5033 fdr < 0.05, 14765 fdr < 0.1. This is great.

### 5. Save differential ATAC results (ASD vs. CTL)
rm(list = ls())
load("1a_03_lme_p_magnitude.rda") # original file name: p_magnitude_lmeDiagnosisASD_ASD_Batch_RegionBA38_RegionBA44.45_tssenrich.score_FRiP.rda
DiffATAC = as.data.frame(cbind(magnitude[,1], p[,1]))
colnames(DiffATAC) = c("ATAC_logFC", "ATAC_p")
DiffATAC$ATAC_FDR = p.adjust(DiffATAC$ATAC_p, method = "fdr")

DiffATAC_location = as.data.frame(str_split_fixed(rownames(DiffATAC), "_", 3))
colnames(DiffATAC_location) = c("chr", "start", "end")
DiffATAC_location$start = as.integer(DiffATAC_location$start)
DiffATAC_location$end = as.integer(DiffATAC_location$end)
DiffATAC = cbind(DiffATAC_location, DiffATAC)

save(DiffATAC, file = "1a_03_DiffATAC.rda") # original file name: DiffATAC.rda


