#!/bin/bash
query_file=$1
file_base=`basename $1 .txt`
output_dir=${2}/${file_base}
mkdir -p ${output_dir}
query_genes=`cat $query_file | perl -pe 'chomp if eof' | tr '\n' '|'`
cat input/meme/all_sc_promoters.fasta | paste -d ':' - - | grep -E $query_genes | tr ':' '\n' > ${output_dir}/${file_base}.fasta
~/meme/bin/ame --verbose 1 --oc ${output_dir}/${file_base} --control input/meme/all_sc_promoters.fasta --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 ${output_dir}/${file_base}.fasta input/meme/JASPAR_CORE_2016_fungi.meme
