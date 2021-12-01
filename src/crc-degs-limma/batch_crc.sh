#!/bin/sh

IN=$1

OUT=$(echo $IN | sed 's@collapsed_norm_data\/mix\/annotation_mix@deg/@g')

echo $IN
echo $OUT

./crc-degs-limma.r ../../data/splitted_samples/stage_low.csv ../../data/splitted_samples/stage_high.csv $IN $OUT
