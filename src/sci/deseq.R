source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')

# Get WTs
dds = dds[,colData(dds)$Strain == 'WT']
# Remove Drug
dds = dds[, colData(dds)$Drug == 'None']
# Remove untrusted
dds = dds[ , colData(dds)$Prep.QC == 'PASSED']

dds = DESeq(dds)

Conditions = unique(colData(dds)$Condition)
reference_condition = which(Conditions == 'None_YPD_None')
Conditions = Conditions[-reference_condition]
res_all = lapply(as.character(Conditions),
                 function(x) results(dds,
						                         contrast=c('Condition', x, 'None_YPD_None')))
