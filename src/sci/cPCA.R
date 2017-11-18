# prepair for cPCA

load('../../input/images/normalized_data.RData')
library(feather)

rownames(t2g) = t2g$target_id
rownames(vlog) = t2g[rownames(vlog),'name']

colnames(vlog) = meta$Experiment
wt = vlog[,meta$Strain == 'WT']
as = vlog[,meta$Strain != 'WT']
write_feather(as.data.frame(wt), 'wt.feather')
write_feather(as.data.frame(as), 'as.feather')
getwd()
