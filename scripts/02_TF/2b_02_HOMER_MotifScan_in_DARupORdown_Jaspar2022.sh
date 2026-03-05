#!/bin/bash

## Goal: Motif enrichment using HOMER

qrsh -l h_data=10G,h_rt=8:00:00
BED_DIR=~/Desktop/ASD_manuscript/14_GitHub/Upload/results/02_TF/ # original directory: /u/project/geschwind/jennybea/ASD_project_2ndbatch/Integrative/9_DiffATACupORdown
cd $BED_DIR
ls *.bed # 2 bed files

PATH=$PATH:/u/project/geschwind/jennybea/Hi-C_tools/HOMER/bin/
KNOWN_motifs=~/Desktop/ASD_manuscript/14_GitHub/Upload/utils/jaspar2022_BrExpTF.motifs # original file: /u/project/geschwind/jennybea/ATAC/NEW_ASD_parietal/Motif/jaspar2022_BrExpTF.motifs

for bed in *.bed
do
  OUT=$BED_DIR/${bed%.bed}_Jaspar2022/ # no need to mkdir, homer will do it.
  if [ ! -d $OUT ]; then
    findMotifsGenome.pl $BED_DIR/$bed hg19 $OUT -size given -nomotif -mknown $KNOWN_motifs & # -nomotif won't find de novo motifs.
  fi
done
# Run for 1hr30min

# Change file names
cd $BED_DIR
for j in *.bed; do
  i=${j%.bed}_Jaspar2022
  mv $BED_DIR/${i}/knownResults.txt $BED_DIR/2a_01_${i}_knownResults.txt
done
