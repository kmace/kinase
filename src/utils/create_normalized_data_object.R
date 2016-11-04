rm(list=ls())
source('../utils/load_libraries.R')
source('../utils/load_functions.R')

t2g = load_transcripts_to_genes()
meta = get_sample_meta()
raw_counts = load_count_matrix(meta,'tophat')

dds = DESeqDataSetFromMatrix(countData = raw_counts, colData = meta, ~1)

colData(dds)$Strain = relevel(factor(colData(dds)$Strain), 'WT')
colData(dds)$Drug = relevel(factor(colData(dds)$Drug), 'None')
colData(dds)$Stress = relevel(factor(colData(dds)$Stress), 'None')
colData(dds)$Media = relevel(factor(colData(dds)$Media), 'YPD')

rld_all <- rlogTransformation(dds, blind=TRUE)
rlog = assays(rld_all)
rlog_meta = colData(rld_all)

save.image('../../input/normalized_data.RData')
