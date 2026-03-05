#!/bin/bash

bam=$1
GENOME_FA=$2
PEAK=$3
blacklist=$4
OUT=$5

TOBIAS ATACorrect --bam $bam --genome $GENOME_FA --peaks $PEAK --blacklist $blacklist --outdir $OUT --cores 8
