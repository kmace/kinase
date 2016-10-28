source('../utils/0.1-data_loading_functions.R')
source('../utils/0.2-data_manipulation_functions.R')
source('../utils/0.3-plotting_functions.R')

t2g = load_transcripts_to_genes()
meta = get_sample_meta()
raw_counts = load_count_matrix(meta,'tophat')

dds_all = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta, ~1)

colData(dds_all)$Strain = relevel(factor(colData(dds_all)$Strain), 'WT')
colData(dds_all)$Drug = relevel(factor(colData(dds_all)$Drug), 'None')
colData(dds_all)$Stress = relevel(factor(colData(dds_all)$Stress), 'None')
colData(dds_all)$Media = relevel(factor(colData(dds_all)$Media), 'YPD')

rld_all <- rlogTransformation(dds_all)

rld_all_filtered = rld_all[,colData(rld_all)$Prep.QC == "PASSED"]

rlogm_all <- assay(rld_all_filtered)

col_lab = factor(paste(colData(rld_all_filtered)$Stress,
                       colData(rld_all_filtered)$Drug,
                       colData(rld_all_filtered)$Media,
                       sep='_' ))
text_lab = factor(paste(colData(rld_all_filtered)$Stress,
                        colData(rld_all_filtered)$Drug,
                        colData(rld_all_filtered)$Media,
                        colData(rld_all_filtered)$Strain,
                       sep='_' ))

#lab = factor(colData(rld_all)$Plate_Code)


colnames(rlogm_all) = text_lab #colData(rld_all)$Stress
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 2
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])


hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
col_lab = col_lab[-leave]
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)

# plot(as.phylo(hc_subset),
#      type = "unrooted",
#      tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))),
#      label.offset=10,
#      cex=.1)


plot(as.phylo(hc_subset),
     type = "fan",
     tip.color = hsv(as.numeric(col_lab)/max(as.numeric(col_lab))),
     label.offset=10,
     cex=.5)

# legend('topright',
#        legend=unique(col_lab),
#        fill=hsv(as.numeric(unique(col_lab))/max(as.numeric(unique(col_lab)))))


plotHclustColors(rlogm_all, colData(rld_all)$Stress)

dds_wt = dds_all[,colData(dds_all)$Strain == 'WT']
dds_wt_noDrug = dds_wt[, colData(dds_wt)$Drug == 'None']

rld <- rlogTransformation(dds_wt_noDrug)
rlogm <- assay(rld)

meanSdPlot(rlogm)
lab = colData(rld)$Stress
plotHclustColors(rlogm, lab)


res_all = lapply(as.character(
            unique(colData(dds_wt_noDrug)$state)
            )[-2], function(x) results(dds_wt_noDrug,
            contrast=c('state', x, 'None_YPD')))
