#!/bin/bash
for folder in /lab/solexa_public/Pincus/*WIGTC*; do 
    file=`basename ${folder}/*annotation.xls`; 
    echo $folder $file >> files.txt; 
    cat ${folder}/$file >> all_meta.txt
done

