source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')
options(mc.cores = 4L)
raw_counts = raw_counts[ ,meta$Strain=='WT']
meta = meta[meta$Strain=='WT',]

meta$sample = meta$Sample_Name
meta$path = meta$kallisto_path
# Get simple:

meta = meta[!(meta$Stress != 'None' & meta$Drug != 'None'),]
meta$Condition= relevel(factor(meta$Condition),"None_YPD_None")

model = model.matrix(~ 1 + Condition, data=meta)
so <- sleuth_prep(meta, model, target_mapping = t2g)
so <- sleuth_fit(so)

for (c in levels(meta$Condition)[-1]){
  so = sleuth_wt(so, which_beta = paste('Condition',c, sep=''))
}
so <- sleuth_wt(so, which_beta = 'ConditionHS')
#so <- sleuth_wt(so, which_beta = 'Drug1NM')
sleuth_live(so)

#results_table <- sleuth_results(so, 'StrainWT')
pinc = sleuth_to_matrix(so,"obs_norm","est_counts")$data
