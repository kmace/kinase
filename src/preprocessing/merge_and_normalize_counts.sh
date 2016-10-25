#!/bin/bash

cd counts

afile=`ls *.txt | head -n 1`

echo "Gene" > geneExpression.xls
cut -f1 $afile >> geneExpression.xls

for a in *.txt; do
    echo $a
    sample=`basename $a _gene_counts.txt`
    echo $sample > temp_sample
    cut -f2 $a >> temp_sample
    paste geneExpression.xls temp_sample > temp_expression;
    mv temp_expression geneExpression.xls;
done
mv geneExpression.xls geneExpression_raw.xls
/nfs/BaRC_Public/BaRC_code/R/DESeq2/DESeq2_normalize_only.R geneExpression_raw.xls geneExpression_all.xls
nf=`cat geneExpression_all.xls | awk '{print NF}' | sort -nu | tail -n 1`
cat geneExpression_all.xls | cut -f 1,$(( $nf - $(( $(( $nf - 1 )) / 2 )) + 1 ))-$nf > geneExpression_normalized.xls
rm temp_sample

Rscript ../../annotate_yeast.R

cd ..