#!/bin/bash
  
parallel -j20 bash run_meme.sh {} ::: ../kinase/output/Strain_Dependent_Genes/*.txt

for a in *_sig*.fasta do 
	b=`basename $a .fasta`; 
	mv $a $b; 
done

for a in *sig; do 
	b=`basename $a _sig`; 
	mkdir $b; 
	mv ${a}* $b; 
done
