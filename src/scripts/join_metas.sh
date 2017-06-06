# Should handle the first line differently,
# also, shold only keep one of the first lines

awk '{print FILENAME, $0}' */meta/*.xls > all_meta.tsv

#for metafile in */meta/*.xls; do 
#    echo $a; 
#    
#done