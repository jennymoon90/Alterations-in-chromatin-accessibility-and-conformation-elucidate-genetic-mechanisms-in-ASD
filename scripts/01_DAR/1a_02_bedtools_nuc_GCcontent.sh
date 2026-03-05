#!/bin/bash
# qrsh -l h_data=10G,h_rt=8:00:00

DIR_result=/Desktop/ASD_manuscript/14_GitHub/Upload/results/01_DAR/
DIR_data=../../data/ # oritinal data name Peaksets_minOverlap5.bed
GENOME=../../utils/ # original directory: /u/project/geschwind/dhg/jennybea/ATAC/UCSC_hg19_rebuildBowtie2myself/
cd $DIR_result

# make bed file for consensus peaks

bedtools nuc -fi $GENOME/hg19.fa -bed $DIR_data/ATAC_peaks.bed > $DIR_result/ATAC_bedtoolsnuc.tsv # original output file name: Peaksets_minOverlap5_bedtoolsnuc.tsv
# hg19.fa.fai not found, generating...
# very fast

vim $DIR_result/ATAC_bedtoolsnuc.tsv
# i -> delete the # in the first row -> esc -> :wq
