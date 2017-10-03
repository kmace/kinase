#!/bin/bash
# This is a script to go from a standard rna-seq experiment all the way through

# One must pass in the location of the solexa_public directory that has all the data

mkdir project_dir
cd project_dir

cp /lab/solexa_public/Pincus/151222_WIGTC-HISEQ2A_C8DB6ACXX/*laneannotation.xls .
/nfs/BaRC_Public/BaRC_code/Perl/download_fq_from_solexa_public/download_fq_from_solexa_public.pl /lab/solexa_public/Pincus/151222_WIGTC-HISEQ2A_C8DB6ACXX fastq

for i in `/bin/ls *fastq/*.fq.gz`
do
   name=`basename $i|perl -pe 's/\fq\.gz//'`
   bsub "tophat --solexa1.3-quals -p 4 --segment-length 20 -I 2500 -G /nfs/genomes/sgd_2010/gtf/Saccharomyces_cerevisiae.R64-1-1.80_chr.gtf -o ${name}_tophat_out /nfs/genomes/sgd_2010/bowtie/sacCer3 $i"
done

mkdir tophat_out_all
for i in `/bin/ls *_tophat_out/accepted_hits.bam`
do
	name=`dirname $i | perl -pe 's/(.*)\_tophat_out/$1/'`
	mv $i tophat_out_all/${name}.bam
done

mkdir tophat_out_all/other_files
mv *_tophat_out tophat_out_all/other_files

mkdir gene_counts
for i in `/bin/ls tophat_out_all/*.bam`
do
	name=`basename $i |perl -pe 's/\.bam//'`
	bsub "htseq-count -f bam -m intersection-strict --stranded=reverse $i /nfs/genomes/sgd_2010/gtf/Saccharomyces_cerevisiae.R64-1-1.80_chr.gtf > gene_counts/${name}.sgd_2010_gene_counts.txt"
done

cd gene_counts
/bin/ls *gene_counts.txt | perl -pe 's/.sgd_2010_gene_counts.txt//g; s/\n/\t/g' | perl -pe 's/\t$/\n/; s/^/Gene\t/;' >| Samples
paste *gene_counts.txt | cut -f1,2,4,6,8,10,12 | head -n -5 >| Gene_counts_all
cat Samples Gene_counts_all >| Gene_counts_all.txt
# Delete intermediate files
rm -f Samples Gene_counts_all

/nfs/BaRC_Public/BaRC_code/R/DESeq2/DESeq2_normalize_only.R Gene_counts_all.txt Gene_counts_all.normalized_too.txt
