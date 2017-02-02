library(DESeq2)
source('0.2-data_manipulation_functions.R')
load('../../input/images/normalized_data.RData')
vsd = vsd[, colData(vsd)$Drug == 'Cocktail']
dat = assay(vsd)
dat = rename_gene_names(dat,t2g)
colnames(data) = colData(vsd)$Condition
meta = colData(vsd)



write.csv(dat,'../../input/nn/data.csv')
write.csv(meta,'../../input/nn/meta.csv')
