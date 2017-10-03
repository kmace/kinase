#!/bin/bash

for file in ../data/Pincus_data/fastq/*.txt.tar.gz; do
  b=`basename $file`
  echo $b
  ../tools/kallisto_linux-v0.42.4/kallisto quant -i ../data/transcripts.idx -b 500 -o ../data/Pincus_data/results/${b} --single -l 321.05 -s 30 -t 12 $file
done
