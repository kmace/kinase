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


dds = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta, ~ Condition + Strain)
dds = dds[,dds$Experimenter=='Kieran']
dds = dds[,dds$Drug == 'Cocktail']


dds = DESeq(dds, parallel = TRUE)
rlt = rlogTransformation(dds, blind = FALSE) # This takes a long long time
rlog = assay(rlt)
meta = colData(rlt)


# dds <- estimateSizeFactors(dds)
# baseMean <- rowMeans(counts(dds, normalized=TRUE))
# sum(baseMean > 1)
# idx <- sample(which(baseMean > 5), 1000)
# dds.sub <- dds[idx, ]
# dds.sub <- estimateDispersions(dds.sub)
# dispersionFunction(dds) <- dispersionFunction(dds.sub)
# vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

save.image('../../input/images/normalized_data.RData')
