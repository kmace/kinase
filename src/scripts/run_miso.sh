
for a in *.sam; do

sam_to_bam --convert 10t_sequence.txt.gz_Aligned.out.sam test.bam
miso --run indexed/ test.sorted.bam --output-dir my_sample_output/ --read-len 50

summarize_miso --summarize-samples SE/control/ SE/control/
summarize_miso --summarize-samples SE/knockdown/ SE/knockdown/
compare_miso --compare-samples SE/control/ SE/knockdown/ SE/comparisons/
