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

meta %<>% filter(Drug == 'Cocktail')
meta %<>% filter(!(Prep.QC == 'MARGINAL' & Sample_Name == '28r'))

raw_counts = load_count_matrix(meta,'star')

raw_counts = rename_gene_names(raw_counts, t2g)



#meta = meta %>% filter(Strain != 'CDC15')
meta$Strain = droplevels(meta$Strain)
meta$Condition = droplevels(meta$Condition)

library("BiocParallel")
register(MulticoreParam(20))

library(DESeq2)

dds = DESeqDataSetFromMatrix(countData = raw_counts[,meta$Sample_Name], colData = meta, ~1)
dds = estimateSizeFactors(dds)
normalized_counts = counts(dds, normalized=TRUE)

# Full design, include everything:
meta_sub = meta[meta$Strain!='CDC15',]
meta_sub$Strain = droplevels(meta_sub$Strain)
design = model.matrix(~Condition + Strain + Condition:Strain, meta_sub)
all.zero <- apply(design, 2, function(x) all(x==0))
idx <- which(all.zero)
design <- design[,-idx]
DESeq2:::checkFullRank(design)

dds_full = DESeq(dds[,dds$Strain!='CDC15'], full = design, parallel = TRUE)

# Normilized data, taking full design into account
vst_full = varianceStabilizingTransformation(dds_full, blind = FALSE)
vlog_full = assay(vst_full)

# Only model the General effects:
dds_reduced = dds
dds_reduced@design = ~ Condition + Strain
dds_reduced = DESeq(dds_reduced, parallel = TRUE)

# One could also use the reduced model, ~Condition + Strain
# to estimate the dispersion. this is what I've been using this whole time come to think of it.
# Normilized data, taking reduced design into account
vst_reduced = varianceStabilizingTransformation(dds_reduced, blind = FALSE)
vlog_reduced = assay(vst_reduced)

# Only model the Condition effects on wt:
dds_wt = dds[,dds$Strain=='WT']
dds_wt@design = ~Condition
dds_wt = DESeq(dds_wt, parallel = TRUE)

# Normilized data, taking wt dispersion into account
dispersionFunction(dds) <- dispersionFunction(dds_wt)
vst = varianceStabilizingTransformation(dds, blind = FALSE)
vlog = assay(vst)

vst_blind = varianceStabilizingTransformation(dds, blind = TRUE)
vlog_blind = assay(vst_blind)

#dds_lrt = dds[,colData(dds)$Strain != 'CDC15']
#dds_lrt = dds_lrt[,colData(dds)$Condition != 'Menadione']
#colData(dds_lrt)$Strain = droplevels(colData(dds_lrt)$Strain)
#colData(dds_lrt)$Condition = droplevels(colData(dds_lrt)$Condition)
#m1 = model.matrix(~Strain*Condition, colData(dds_lrt))
#m2 = model.matrix(~Strain + Condition, colData(dds_lrt))

#all.zero <- apply(m1, 2, function(x) all(x==0))
#idx <- which(all.zero)
#m1 <- m1[,-idx]
#dds_lrt = DESeq(dds_lrt, test = 'LRT', full = m1, reduced = m2, betaPrior = FALSE)



meta = as_tibble(meta)
t2g$Gene = t2g$target_id
t2g = as_tibble(t2g)


save(file = '../../intermediate/images/normalized_data.RData',
     list = c('meta', 't2g', 'raw_counts', 'normalized_counts',
              'dds', 'vst', 'vlog', 'vst_blind', 'vlog_blind', # dispersions from WT
              'dds_full', 'vst_full', 'vlog_full',
              'dds_reduced', 'vst_reduced', 'vlog_reduced',
              'dds_wt',
              'dds', 'vst_blind', 'vlog_blind'
     ))
