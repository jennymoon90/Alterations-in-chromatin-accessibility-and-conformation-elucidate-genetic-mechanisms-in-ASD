qrsh -l h_data=10G,h_rt=12:00:00

### Note: In the following script, I combined all ASD samples into one .bw and all CTL samples into one .bw file. The TF_bound.bed files are also with the combined .bw file. But to validate the footprint shapes, one may run with any one sample for the TF (replace BINDetect_Jaspar2022/ with indBINDetect_Jaspar2022/)

DIR=~/Desktop/ASD_manuscript/14_GitHub/Upload/ # original: /u/project/geschwind/jennybea/ASD_project_2ndbatch/ATAC/
OUT=results/02_TF/2c_01_TOBIAS/ #original: $DIR/10_TOBIAS/
OUT_plot=plots/02_TF/
Debug=$DIR/debug/
SCRIPT=scripts/02_TF/ # original: $DIR/scripts/
PEAK=$DIR/data/ATAC_peaks.bed # original: $DIR/4_ConsensusPeaks/Diffbind/Peaksets_minOverlap5.bed
GENOME_FA=$DIR/utils/hg19.fa # original: /u/project/geschwind/dhg/jennybea/ATAC/UCSC_hg19_rebuildBowtie2myself/hg19.fa
MOTIF=$DIR/utils/Jaspar2022motifs.jaspar # original: /u/project/geschwind/jennybea/Resource/TFmotifs/Jaspar2022motifs.jaspar.

## step3. BINDetect
mkdir -p $OUT/BINDetect_Jaspar2022/
cd $OUT/BINDetect_Jaspar2022/
qsub -V -cwd -o $Debug -e $Debug -N BINDetect_ASD \
-l h_data=15G,h_rt=8:00:00,highp -pe shared 8 \
$SCRIPT/2c_04_agg_BINDetect.sh $MOTIF \
$OUT/ScoreBigwig/CTL_footprints.bw $OUT/ScoreBigwig/ASD_footprints.bw \
$GENOME_FA $PEAK $OUT/BINDetect_Jaspar2022/ CTL ASD # original script: 10_2_3_p_BINDetect.sh
# Run for 20min.

## step4. PlotAggregate
mkdir -p $OUT_plot/2c_04_PlotAggregate/ # original: $OUT/PlotAggregate_Jaspar2022/
cd $OUT/BINDetect_Jaspar2022/

while read -r line;
do
  TFname="$line"
  TFdir=${TFname}_${TFname}
  tmp=$(echo $TFname | cut -d"-" -f3) # should be TFshort=

  if [[ ! -f $OUT_plot/2c_04_PlotAggregate/${TFname}_footprint_comparison_bound.pdf ]]; then
    echo $TFname
    tfbs_asd=$OUT/BINDetect_Jaspar2022/${TFdir}/beds/${TFdir}_ASD_bound.bed
    tfbs_ctl=$OUT/BINDetect_Jaspar2022/${TFdir}/beds/${TFdir}_CTL_bound.bed

    TOBIAS PlotAggregate --TFBS $tfbs_ctl $tfbs_asd --signals $OUT/ATACorrect/CTL_merged_downsizedto100M_corrected.bw $OUT/ATACorrect/ASD_merged_downsizedto100M_corrected.bw --output $OUT_plot/2c_04_PlotAggregate/${TFname}_footprint_comparison_bound.pdf --share_y both --plot_boundaries  --signal-labels CTL ASD --TFBS-labels ${TFshort}_CTL_bound ${TFshort}_ASD_bound
  fi
done < $OUT/../2c_03_PlotAggregate.txt # original: $OUT/PlotAggregate_Jaspar2022/DiffTFbindingscoreFromlmeVolcano_Jaspar2022TFmotif_for_PlotAggregate.txt
# Run for 15min.
