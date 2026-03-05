#!/bin/bash
qrsh -l h_data=10G,h_rt=8:00:00
TOBIAS
# Refer to https://github.com/loosolab/TOBIAS for pipeline

## step1. Bias correction of ATAC-seq reads
DIR=~/Desktop/ASD_manuscript/14_GitHub/Upload/ # original: /u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/
BAM=$DIR/data/ATAC_bam_UniqueAligned/ # original: $DIR/2_MappingToPeakcalling/bam_UniqueAligned/
PEAK=$DIR/data/ATAC_peaks.bed # original: $DIR/4_ConsensusPeaks/Diffbind/Peaksets_minOverlap5.bed
GENOME_FA=$DIR/utils/hg19.fa # original: /u/project/geschwind/dhg/jennybea/ATAC/UCSC_hg19_rebuildBowtie2myself/hg19.fa
blacklist=$DIR/utils/hg19.blacklist.bed # original /u/project/geschwind/jennybea/ATAC/blacklist/hg19.blacklist.bed # hg38.blacklist.bed
OUT=results/02_TF/2c_01_TOBIAS/ #original: $DIR/10_TOBIAS/
mkdir -p $OUT/ind_ATACorrect/
Debug=$DIR/debug/
SCRIPT=scripts/02_TF/ # original: $DIR/scripts/

ASD_good_samples=("2GGP-1"  "2NA6-1"  "35BV-1"  "75QW-1"  "B1793-1" "B1860"   "B4334-1" "B4334-2" "B4337_4" "B4721A"  "B4721B"  "B4787"   "B5000B"  "B5000C"  "B5144"   "B5163-2" "B5242"   "B5297-1" "B5309"   "B5352_1" "B5558-2" "B5569-1" "B5666"   "B5718B"  "B5718D"  "B5813A"  "B5813B" "B6200"   "CQ56-1"  "CQ56-2"  "EUFW-1"  "FHN9-2"  "G7YU"    "JU7U-1"  "KMC3-1"  "M9H3-1" "PH2A-1"  "VPSP")

cd $BAM
for s in "${ASD_good_samples[@]}"; do
  if [[ ! -f $OUT/ind_ATACorrect/${s}_unique_corrected.bw ]]; then
    qsub -V -cwd -o $Debug -e $Debug -N ATACorrect_$s \
    -l h_data=4G,h_rt=4:00:00,highp -pe shared 8 \
    $SCRIPT/2c_01_runTOBIAS_st1_ATACorrect.sh $BAM/${s}_unique.bam $GENOME_FA $PEAK $blacklist $OUT/ind_ATACorrect/ # original script 10_2_1_p_ATACorrect.sh
  fi
done
# Run for 7hrs

## step2. ScoreBigwig
mkdir -p $OUT/ind_ScoreBigwig/

# To view TF binding change with age and diagnosis, calculate the binding score and plot. Run the following codes:
cd $OUT/ind_ATACorrect/
for i in *_unique_corrected.bw; do
  s=${i%_unique_corrected.bw}
  if [[ ! -f $OUT/ind_ScoreBigwig/${s}_footprints.bw ]]; then
    qsub -hold_jid ATACorrect_$s -V -cwd -o $Debug -e $Debug -N ScoreBigwig_$s \
    -l h_data=15G,h_rt=4:00:00,highp -pe shared 8 \
    $SCRIPT/2c_01_runTOBIAS_st2_ScoreBigwig.sh $OUT/ind_ATACorrect/${i} $PEAK $OUT/ind_ScoreBigwig/${s}_footprints.bw # original script 10_2_2_p_ScoreBigwig.sh
  fi
done
# Run for 5.5hrs

## step3. BINDetect on single condition (single sample)
MOTIF=$DIR/utils/Jaspar2022motifs.jaspar # original: /u/project/geschwind/jennybea/Resource/TFmotifs/Jaspar2022motifs.jaspar.

cd $OUT/ind_ATACorrect/
for i in *_unique_corrected.bw; do
  s=${i%_unique_corrected.bw}
  mkdir -p $OUT/ind_BINDetect_Jaspar2022/${s}/
  if [[ ! -f $OUT/ind_BINDetect_Jaspar2022/${s}/bindetect_results.txt ]]; then
    qsub -hold_jid ScoreBigwig_$s -V -cwd -o $Debug -e $Debug -N BINDetect_$s \
    -l h_data=10G,h_rt=8:00:00,highp -pe shared 8 \
    $SCRIPT/2c_01_runTOBIAS_st3_BINDetect_ind.sh $MOTIF \
    $OUT/ind_ScoreBigwig/${s}_footprints.bw \
    $GENOME_FA $PEAK $OUT/ind_BINDetect_Jaspar2022/${s}/ $s # original script 10_3_3_p_BINDetect_ind.sh
  fi
done
# Run for 1hr10min

# to plot and see whether there is any interesting trajectory changes b/w ASD and CTL
cd $OUT/ind_BINDetect_Jaspar2022/
mkdir -p all_bindetect_results/
for s in *; do
  if [[ $s != "all_bindetect_results" && ! -f $OUT/ind_BINDetect_Jaspar2022/all_bindetect_results/${s}_bindetect_results.txt ]]; then
    ln -s $OUT/ind_BINDetect_Jaspar2022/$s/bindetect_results.txt $OUT/ind_BINDetect_Jaspar2022/all_bindetect_results/${s}_bindetect_results.txt
  fi
done

cd all_bindetect_results/; ls
