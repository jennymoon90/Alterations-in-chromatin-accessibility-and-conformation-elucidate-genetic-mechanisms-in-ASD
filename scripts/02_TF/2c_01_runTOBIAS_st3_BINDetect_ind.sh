#!/bin/bash

motif=$1
bw=$2
GENOME_FA=$3
PEAK=$4
out=$5
name=$6

TOBIAS BINDetect --motifs $motif --signals $bw --genome $GENOME_FA --peaks $PEAK \
--outdir $out --cond_names $name --cores 8
