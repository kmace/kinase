#!/bin/bash
query_file=$1
file_base=`basename $1 .txt`
output_dir=${2}

mkdir -p ${output_dir}/fasta
mkdir -p ${output_dir}/meme
mkdir -p ${output_dir}/tfs

query_genes=`cat $query_file | perl -pe 'chomp if eof' | tr '\n' '|'`
cat input/meme/all_sc_promoters.fasta | paste -d ':' - - | grep -E $query_genes | tr ':' '\n' > ${output_dir}/fasta/${file_base}.fasta

~/meme/bin/ame --verbose 1 --oc ${output_dir}/meme/${file_base} --control input/meme/all_sc_promoters.fasta --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 ${output_dir}/fasta/${file_base}.fasta input/meme/JASPAR_CORE_2016_fungi.meme

src/scripts/meme/get_TFs.sh ${output_dir}/meme/${file_base}/ame.html > ${output_dir}/tfs/${file_base}_tfs.txt
