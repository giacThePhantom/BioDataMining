#!/bin/sh


IN=$1

OUT=$(echo $IN | sed 's/norm_data/collapsed_norm_data/g')

echo $IN
echo $OUT

./collapse_genes.r $IN $OUT
