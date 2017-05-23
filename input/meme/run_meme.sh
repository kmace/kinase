#!/bin/bash
query_file=$1
file_base=`basename $1 .txt`
query_genes=`cat $query_file | perl -pe 'chomp if eof' | tr '\n' '|'`
cat all_sc_promoters.fasta | paste -d ':' - - | grep -E $query_genes | tr ':' '\n' > ${file_base}.fasta
../meme/bin/ame --verbose 1 --oc ./${file_base} --control all_sc_promoters.fasta --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 ${file_base}.fasta JASPAR_CORE_2016_fungi.meme
