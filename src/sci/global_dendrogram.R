source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')

passed = rlog_meta$Prep.QC == "PASSED"
rlog = rlog[ ,passed]
rlog_meta = rlog_meta[passed, ]
rlog <- assay(rld_all_filtered)

col_lab = factor(paste(rlog_meta$Stress,
											 rlog_meta$Drug,
											 rlog_meta$Media,
											 sep='_' ))

text_lab = factor(paste(rlog_meta$Stress,
												rlog_meta$Drug,
												rlog_meta$Media,
												rlog_meta$Strain,
											 sep='_' ))

colnames(rlog) = text_lab #colData(rld_all)$Stress
dst <- dist(t(rlog))
hc = hclust(dst)

# use this code to remove outliers:
# ave_dist = apply(as.matrix(dist_all),2,mean)
# num_to_remove = 2
# leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
# hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
# col_lab = col_lab[-leave]




plot(as.phylo(hc_subset),
		 type = "fan",
		 tip.color = hsv(as.numeric(col_lab)/max(as.numeric(col_lab))),
		 label.offset=10,
		 cex=.5)

# legend('topright',
#				legend=unique(col_lab),
#				fill=hsv(as.numeric(unique(col_lab))/max(as.numeric(unique(col_lab)))))



meanSdPlot(rlog)

lab = colData(rld)$Stress
plotHclustColors(rlog, lab)
