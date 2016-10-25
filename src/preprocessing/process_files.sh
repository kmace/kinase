#!/bin/bash                                                                                                                                         # Get the location of the files and desired project name 
source_location=${1%/}
project_name=${2%/}
mkdir $project_name
cd $project_name
# Make holding folders
mkdir fastq tophat_output bam counts meta
# Copy metadata across
cp ${source_location}/*laneannotation.xls meta/
Rscript ../parseMeta.R meta/*laneannotation.xls

# Loop through meta data
jobnum=0
while read STR; do
    sample=$(echo $STR | cut -f1 -d,)
    file=$(echo $STR | cut -f2 -d,)
    echo sending $file for processing
    bsub -J "${project_name}_lib${jobnum}" "../run_one_library.sh $source_location $file $sample"
    jobnum=$((jobnum + 1))
done <meta/meta.tsv


bsub -J "merge_libs" -w "done(${project_name}_lib*)" "../merge_and_normalize_counts.sh"
cd ..



