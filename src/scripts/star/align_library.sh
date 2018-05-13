#!/bin/bash
for i in $(find ../wi.mit.edu/161202_WIGTC-HISEQB_CA4URACXX/fastq/* -name "*.txt.gz"); do
    echo $i;
    hand=`basename $i`
    STAR --genomeDir  genomeDir --readFilesIn $i --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts  --runThreadN 10 --readFilesCommand zcat --outFileNamePrefix ./star_output/${hand}_;
done
