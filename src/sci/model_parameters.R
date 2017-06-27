rm(list=ls())

load('../../input/images/normalized_data.RData')

library(tidyr)
library(dplyr)
library(broom)
library(tibble)

vlog = vlog - rowMeans(vlog)
row_sd = apply(vlog,1,sd)
vlog = apply(vlog,2,function(x) x / row_sd)

vlog = as.data.frame(vlog)
vlog$Gene = rownames(vlog)
measurements = gather(vlog, Experiment, Expression, -Gene)
# should be a left join on the meta subsetted
measurements = measurements %>% 
  mutate(Strain = meta[Experiment, ]$Strain,
         Condition = meta[Experiment, ]$Condition) %>%
  select(-Experiment)
              
fits = measurements %>% 
  dplyr::group_by(Gene) %>%
  dplyr::do(fitGene = lm(Expression ~ Strain + Condition, data = .)) %>%
  broom::tidy(fitGene)

parameters = fits %>% select(Gene, term, estimate) %>% spread(term, estimate) %>% ungroup()
rownames(parameters) = parameters$Gene
parameters = parameters[,-Gene]

parameter_matrix = as.matrix(select(parameters, -Gene))
col_meam = apply(mat,2,mean)
#m = apply(mat,1,function(x) x - col_m)
#m = t(m)
#heatmap3(m[m[,10]>1,-10], labRow = NA, scale='none')
heatmap3(mat, labRow = NA, scale='none')

