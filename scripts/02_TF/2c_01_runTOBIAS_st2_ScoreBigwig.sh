#!/bin/bash

bw=$1
PEAK=$2
out=$3

TOBIAS FootprintScores --signal $bw --regions $PEAK --output $out --cores 8
