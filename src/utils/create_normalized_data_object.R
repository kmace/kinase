rm(list=ls())
source('../utils/load_libraries.R')
source('../utils/load_functions.R')

t2g = load_transcripts_to_genes()
meta = get_sample_meta()

meta$Strain = relevel(factor(meta$Strain), 'WT')
meta$Drug = relevel(factor(meta$Drug), 'None')
meta$Stress = relevel(factor(meta$Stress), 'None')
meta$Media = relevel(factor(meta$Media), 'YPD')
meta = meta %>% mutate(Condition = if_else(Stress=='None',
                                           as.character(Media),
                                           as.character(Stress)
                                           ))
meta$Condition = relevel(factor(meta$Condition), 'YPD')

raw_counts = load_count_matrix(meta,'star')

subset = meta$Experimenter == 'Kieran' & meta$Drug == 'Cocktail'
meta = meta[subset, ]
raw_counts = raw_counts[, subset]
rm(subset) # clean unused variables

dds = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta, ~ Condition + Strain)


dds = DESeq(dds, parallel = TRUE)

# This takes an incredibly long time
# rlt = rlogTransformation(dds, blind = FALSE)
# rlog = assay(rlt)
# meta = colData(rlt)

# This takes a very ling time
vst = varianceStabilizingTransformation(dds, blind = FALSE)
vlog = assay(vst)
meta = colData(vst)

# This is pretty quick
# baseMean <- rowMeans(counts(dds, normalized=TRUE))
# idx <- sample(which(baseMean > 5), 1000)
# dds.sub <- dds[idx, ]
# dds.sub <- estimateDispersions(dds.sub)
# dispersionFunction(dds) <- dispersionFunction(dds.sub)
# rm(dds.sub)
# vst <- varianceStabilizingTransformation(dds, blind=FALSE)
# vlog = assay(vst)
# meta = colData(vst)

save.image('../../input/images/normalized_data.RData')
