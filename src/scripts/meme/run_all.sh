#!/bin/bash

src_dir=$1
output_dir=$2

mkdir -p ${output_dir}

parallel -j20 bash src/scripts/meme/run_meme.sh {} ${output_dir} ::: ${src_dir}/*.txt
