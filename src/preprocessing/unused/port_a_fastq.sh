#!/bin/bash
#bsub "../scripts/port_a_fastq.sh $source_location $file $sample"
source_location=$1
source_filename=$2
untarred_filename=`basename ${source_filename} .tar.gz`
sample_name=$3

cd fastq
tar -xf ${source_location}/QualityScore/${source_filename} # should go to local -C . fastq/${sample_name}_sequence.txt.tar.gz
mv ./lab/solexa_public/Pincus/*/QualityScore/${untarred_filename} ${sample_name}_sequence.txt
gzip ${sample_name}_sequence.txt
cd ..
