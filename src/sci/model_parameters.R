rm(list=ls())

load('../../input/images/normalized_data.RData')

library(tidyr)
library(purrr)
library(dplyr)
library(broom)
library(tibble)
library(heatmap3)
library(ggplot2)
#library(multidplyr)

#vlog = vlog - rowMeans(vlog)
#row_sd = apply(vlog,1,sd)
#vlog = apply(vlog,2,function(x) x / row_sd)

vlog = as.data.frame(vlog)
vlog$Gene = rownames(vlog)
measurements = gather(vlog, Sample_Name, Expression, -Gene)
# should be a left join on the meta subsetted
meta_sub = as.data.frame(meta) %>% select(Sample_Name, Condition, Strain)

measurements = left_join(measurements, meta_sub)

lin_mod = function(formula) {
  function(data,...){
  map(data,~lm(formula, data = .))
  }
}

list_model <- list(Strain_model= Expression ~ Strain,
                   Condition_model= Expression ~ Condition,
                   Full_model= Expression ~ Strain + Condition) %>%
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
  mutate_at(.vars=vars(ends_with('aug')), .funs = funs(resid = map(.,'.resid'))) %>%
  arrange(Gene)
save.image('../../input/images/model_parameters.RData')
