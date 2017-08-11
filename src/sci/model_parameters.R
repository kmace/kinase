rm(list=ls())

load('../../input/images/normalized_data.RData')

library(tidyr)
library(purrr)
library(dplyr)
library(broom)
library(tibble)
library(heatmap3)
library(ggplot2)
#library(magrittr)
#library(multidplyr)

#vlog = vlog - rowMeans(vlog)
#row_sd = apply(vlog,1,sd)
#vlog = apply(vlog,2,function(x) x / row_sd)

vlog = as.data.frame(vlog)
vlog$Gene = rownames(vlog)
measurements = gather(vlog, Sample_Name, Expression, -Gene)
#measurements %<>% as_tibble()

# should be a left join on the meta subsetted
meta_sub = as.data.frame(meta) %>% select(Sample_Name, Condition, Strain)

measurements = left_join(measurements, meta_sub)

lin_mod = function(formula) {
  function(data, ...){
    map(data, ~ lm(formula, data = .))
  }
}

list_model <- list(Strain_model = Expression ~ Strain,
                   Condition_model = Expression ~ Condition,
                   Full_model = Expression ~ Strain + Condition,
                   Base_model = Expression ~ 1) %>%
  lapply(lin_mod)

genes = measurements %>%
  group_by(Gene) %>%
  mutate(Gene_range = diff(range(Expression)),
         Gene_sd = sd(Expression),
         Gene_mean = mean(Expression),
         Norm_Expression = (Expression - Gene_mean) / Gene_sd) %>%
  filter(Gene_range > 1) %>%
  nest() %>%
  mutate_at(.vars=("data"),.funs=list_model) %>%
  mutate_at(.vars=vars(ends_with('model')), .funs = funs(aug = map(.,augment))) %>%
  mutate_at(.vars=vars(ends_with('aug')), .funs = funs(resid = map(.,'.std.resid'))) %>%
  arrange(Gene)


res_matrix = genes %>%
  select(Gene, data, ends_with('resid')) %>%
  unnest() %>%
  group_by(Strain, Gene, Condition) %>%
  # should only apply to WT
  summarize_at(.vars=vars(ends_with('resid')), mean) %>%
  ungroup() %>%
  mutate(sample_id = paste(Condition,Strain, sep = '_')) %>%
  select(Gene, sample_id, ends_with('resid')) %>%
  rename(Strain = Strain_model_aug_resid,
         Condition = Condition_model_aug_resid,
         Full = Full_model_aug_resid,
         Base = Base_model_aug_resid) %>%
  gather(residual_type, residual, -sample_id, -Gene)
#%>% select(residual, Gene, meta) %>% spread(meta, residual)

get_matrix = function(data, type='Strain'){
  cm = data %>%
    filter(residual_type == type) %>%
    select(-residual_type) %>%
    spread(key = sample_id, value = residual, drop=TRUE)
  rownames(cm) = cm$Gene
  cm  = as.matrix(select(cm, -Gene))
  return(cm)
}

base = get_matrix(res_matrix, type='Base')
strain = get_matrix(res_matrix, type='Strain')
condition = get_matrix(res_matrix, type='Condition')
full = get_matrix(res_matrix, type='Full')


colnames(t2g)[1] = 'Gene'

genes = left_join(genes, t2g)
save.image('../../input/images/model_parameters.RData')
