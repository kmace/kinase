#!/bin/bash
source_location=$1
source_filename=$2
untarred_filename=`basename ${source_filename} .tar.gz`
sample_name=$3

#updated to extract the file out of tar ball
cd fastq
tar -xf ${source_location}/QualityScore/${source_filename} # should go to local -C . fastq/${sample_name}_sequence.txt.tar.gz
mv ./lab/solexa_public/Pincus/*/QualityScore/${untarred_filename} ${sample_name}_sequence.txt
gzip ${sample_name}_sequence.txt
cd ..
tophat --solexa1.3-quals -p 4 --segment-length 20 -I 2500 -G /nfs/genomes/sgd_2010/gtf/Saccharomyces_cerevisiae.R64-1-1.80_chr.gtf -o tophat_output/${sample_name} /nfs/genomes/sgd_2010/bowtie/sacCer3 fastq/${sample_name}_sequence.txt.tar.gz
mv tophat_output/${sample_name}/accepted_hits.bam bam/${sample_name}.bam
htseq-count -f bam -m intersection-strict --stranded=reverse bam/${sample_name}.bam /nfs/genomes/sgd_2010/gtf/Saccharomyces_cerevisiae.R64-1-1.80_chr.gtf > counts/${sample_name}_gene_counts.txt
