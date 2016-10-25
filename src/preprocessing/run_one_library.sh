#!/bin/bash
source_location=$1
source_filename=$2
sample_name=$3

cp ${source_location}/QualityScore/${source_filename} fastq/${sample_name}_sequence.txt.tar.gz
tophat --solexa1.3-quals -p 4 --segment-length 20 -I 2500 -G /nfs/genomes/sgd_2010/gtf/Saccharomyces_cerevisiae.R64-1-1.80_chr.gtf -o tophat_output/${sample_name} /nfs/genomes/sgd_2010/bowtie/sacCer3 fastq/${sample_name}_sequence.txt.tar.gz
mv tophat_output/${sample_name}/accepted_hits.bam bam/${sample_name}.bam
htseq-count -f bam -m intersection-strict --stranded=reverse bam/${sample_name}.bam /nfs/genomes/sgd_2010/gtf/Saccharomyces_cerevisiae.R64-1-1.80_chr.gtf > counts/${sample_name}_gene_counts.txt
