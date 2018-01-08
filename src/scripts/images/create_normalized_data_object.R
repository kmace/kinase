rm(list=ls())
source('../utils/load_libraries.R')
source('../utils/load_functions.R')
library(magrittr)

t2g = readr::read_csv('../../intermediate/t2g.csv')#load_transcripts_to_genes()
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

throw = read_csv('../../intermediate/Samples_to_thow_out.csv') %>% left_join(meta)

meta %<>% filter(!(Sample_Name %in% throw$Sample_Name))


raw_counts = load_count_matrix(meta,'star')

raw_counts = rename_gene_names(raw_counts, t2g)



#meta = meta %>% filter(Strain != 'CDC15')
meta$Strain = droplevels(meta$Strain)
meta$Condition = droplevels(meta$Condition)

library("BiocParallel")
register(MulticoreParam(20))

library(DESeq2)

#se = SummarizedExperiment(assays=list(counts=as.matrix(raw_counts)),
#                          colData=meta)

dds = DESeqDataSetFromMatrix(countData = raw_counts[,meta$Sample_Name], colData = meta, ~1)
dds = estimateSizeFactors(dds)
normalized_counts = counts(dds, normalized=TRUE)

# Full design, include everything except CDC15 which is missing its YPD condition:
dds_full = dds
dds_full = dds_full[,dds_full$Strain!='CDC15']
colData(dds_full)$Strain = droplevels(colData(dds_full)$Strain)
full_design = model.matrix(~Condition + Strain + Condition:Strain, colData(dds_full))
all.zero <- apply(full_design, 2, function(x) all(x==0))
idx <- which(all.zero)
full_design = full_design[,-idx]
rm(all.zero, idx)
DESeq2:::checkFullRank(full_design)

dds_full = DESeq(dds_full, full = full_design, parallel = TRUE)
vst_full = varianceStabilizingTransformation(dds_full, blind = FALSE)
vlog_full = assay(vst_full)

# Model Conditions, and then remainder as strain:
dds_conditionFirst = dds
conditionFirst_design = model.matrix(~Condition+Condition:Strain, meta)
all.zero <- apply(conditionFirst_design, 2, function(x) all(x==0))
idx <- which(all.zero)
conditionFirst_design = conditionFirst_design[,-idx]
DESeq2:::checkFullRank(conditionFirst_design)

dds_conditionFirst = DESeq(dds_conditionFirst, full = conditionFirst_design, parallel = TRUE)
vst_conditionFirst = varianceStabilizingTransformation(dds_conditionFirst, blind = FALSE)
vlog_conditionFirst = assay(vst_conditionFirst)

# Model Strain, and then remainder as condition:
dds_strainFirst = dds
dds_strainFirst = dds_strainFirst[,dds_strainFirst$Strain!='CDC15']
colData(dds_strainFirst)$Strain = droplevels(colData(dds_strainFirst)$Strain)
strainFirst_design = model.matrix(~Strain+Strain:Condition, colData(dds_strainFirst))
all.zero <- apply(strainFirst_design, 2, function(x) all(x==0))
idx <- which(all.zero)
strainFirst_design = strainFirst_design[,-idx]
rm(all.zero, idx)
DESeq2:::checkFullRank(strainFirst_design)

dds_strainFirst= DESeq(dds_strainFirst, full = strainFirst_design, parallel = TRUE)
vst_strainFirst = varianceStabilizingTransformation(dds_strainFirst, blind = FALSE)
vlog_strainFirst = assay(vst_strainFirst)

# Only model the General effects:
dds_reduced = dds
reduced_design = model.matrix(~ Condition + Strain, colData(dds_reduced))
dds_reduced@design = ~ Condition + Strain
dds_reduced = DESeq(dds_reduced, full = reduced_design, parallel = TRUE)

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

# LRT test between full and reduced.

dds_lrt = dds
dds_lrt = dds_lrt[,colData(dds_lrt)$Strain != 'CDC15']
#dds_lrt = dds_lrt[,colData(dds)$Condition != 'Menadione']
colData(dds_lrt)$Strain = droplevels(colData(dds_lrt)$Strain)
colData(dds_lrt)$Condition = droplevels(colData(dds_lrt)$Condition)

dds_lrt = DESeq(dds_lrt, test = 'LRT',
                full = full_design,
                reduced = reduced_design[meta$Strain!='CDC15', -which(colnames(reduced_design) == "StrainCDC15")],
                parallel = TRUE)



meta = as_tibble(meta)
t2g$Gene = t2g$target_id
t2g = as_tibble(t2g)


save(file = '../../intermediate/images/normalized_data.RData',
     list = c('meta', 't2g', # Meta data
              'raw_counts', 'normalized_counts', # Count data
              'dds_wt', # only WT data
              'dds', 'vst', 'vlog', 'vst_blind', 'vlog_blind', # dispersions for dds from WT only (No tests in dds)
              'dds_reduced', 'vst_reduced', 'vlog_reduced', 'reduced_design', # Simple model, no interactions
              'dds_conditionFirst', 'vst_conditionFirst', 'vlog_conditionFirst', 'conditionFirst_design',# Model YPD -> condition, then the strains
              'dds_strainFirst', 'vst_strainFirst', 'vlog_strainFirst', 'strainFirst_design',# Model WT -> strain, then the conditions
              'dds_full', 'vst_full', 'vlog_full', 'full_design',# Model strain and condition as parallel processes, and then the interaction as additional signal after the two parallel processes have happened.
              'dds_lrt' # LRT to test a full model vs the reduced model (is there at least one significant interaction term for a given gene.)
     ))
