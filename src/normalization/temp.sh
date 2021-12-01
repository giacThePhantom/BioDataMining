#!/bin/sh

ls ../../data/norm_data/mix | xargs -P 4 -I {} ./collapse_batch.sh ../../data/norm_data/mix/{}

./collapse_batch.sh ../../data/norm_data/annotation_os_frma.txt
