#!/bin/bash                                                                                                                                         
# Get the location of the files and desired project name 
source_location=${1%/}
project_name=${2%/}

cd ${project_name}

# Loop through meta data

while read STR; do
    sample=$(echo $STR | cut -f1 -d,)
    file=$(echo $STR | cut -f2 -d,)
    echo sending $file for processing
    bsub "../scripts/port_a_fastq.sh $source_location $file $sample"
done <meta/meta.tsv



cd ..



