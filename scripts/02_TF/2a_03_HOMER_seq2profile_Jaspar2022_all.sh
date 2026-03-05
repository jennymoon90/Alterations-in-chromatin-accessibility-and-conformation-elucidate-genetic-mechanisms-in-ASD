#!/bin/bash

qrsh -l h_data=10G,h_rt=8:00:00

PATH=$PATH:/u/project/geschwind/jennybea/Hi-C_tools/HOMER/bin/
cd ~/Desktop/ASD_manuscript/14_GitHub/Upload/scripts/02_TF/ # original: /u/project/geschwind/jennybea/ATAC/NEW_ASD_parietal/Motif/
bash 2a_02_HOMER_seq2profile_Jaspar2022.sh # original script: HOMER_seq2profile_Jaspar2022.sh

mkdir -p ../../utils/res_seq2profile_Jaspar2022
mv *.motif ../../utils/res_seq2profile_Jaspar2022/
cd ../../utils/res_seq2profile_Jaspar2022
cat *.motif > Jaspar2022_BrExp_withHOMERlogodds.motif
