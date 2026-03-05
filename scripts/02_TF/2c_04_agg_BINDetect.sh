#!/bin/bash

motif=$1
bw1=$2
bw2=$3
GENOME_FA=$4
PEAK=$5
out=$6
name1=$7
name2=$8

TOBIAS BINDetect --motifs $motif --signals $bw1 $bw2 --genome $GENOME_FA --peaks $PEAK \
--outdir $out --cond_names $name1 $name2 --cores 8
