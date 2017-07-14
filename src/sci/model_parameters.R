rm(list=ls())

load('../../input/images/normalized_data.RData')

library(tidyr)
library(dplyr)
library(broom)
library(tibble)
library(heatmap3)
library(ggplot2)

#vlog = vlog - rowMeans(vlog)
#row_sd = apply(vlog,1,sd)
#vlog = apply(vlog,2,function(x) x / row_sd)

vlog = as.data.frame(vlog)
vlog$Gene = rownames(vlog)
measurements = gather(vlog, Sample_Name, Expression, -Gene)
# should be a left join on the meta subsetted
meta_sub = as.data.frame(meta) %>% select(Sample_Name, Condition, Strain)

measurements = left_join(measurements, meta_sub)

measurements = measurements %>%
  group_by(Gene) %>%
  mutate(Gene_range = diff(range(Expression)),
         Gene_sd = sd(Expression),
         Gene_mean = mean(Expression),
         Norm_Expression = (Expression - Gene_mean) / Gene_sd) %>%
  ungroup() %>%
  arrange(Gene)

measurements = measurements %>% filter(Gene_range > 1)

fits = measurements %>%
  dplyr::group_by(Gene) %>%
  dplyr::do(fullFit = lm(Expression ~ Strain + Condition, data = .),
            strainFit = lm(Expression ~ Strain, data = .),
            conditionFit = lm(Expression ~ Condition, data = .),
            basicFit = lm(Expression ~ 1, data = .)) %>%
  ungroup() %>%
  tidyr::gather(Model_Type, Model, -Gene) %>%
  rowwise() 

condition_residuals = measurements %>%
  dplyr::group_by(Gene) %>%
  dplyr::do(conditionFit = lm(Expression ~ Condition, data = .)) %>%
  ungroup() %>%
  tidyr



parameters = fits %>%
  broom::tidy(Model)

strength = fits %>%
  broom::glance(Model)

strength %>% 
  filter(Model_Type != 'basicFit') %>% 
  ggplot(aes(x = r.squared, color=Model_Type, fill = Model_Type)) + 
  geom_density(alpha=.5)

sq = strength %>% ungroup() %>% select(name, Model_Type, r.squared) %>%
  filter(Model_Type != 'basicFit') %>% 
  spread(Model_Type, r.squared)

all_models = strength %>%
             filter(Model_Type != 'basicFit') %>%
             ggplot(aes(x = r.squared, 
                        fill = Model_Type)) + 
             geom_density(alpha = 0.5) + 
             ggtitle(label = 'Comparing fits of Models', subtitle = 'Percentage of Variance Explained by each model')

all_models

strain_vs_condition = ggplot(sq, 
                             aes(x = conditionFit, 
                                 y = strainFit, 
                                 name = name)) + 
                    geom_point() + 
                    ggtitle(label = 'Comparing fits of Partial Models', subtitle = 'Percentage of Variance Explained by each model')

additive_vs_full =  ggplot(sq, 
                           aes(x = conditionFit + strainFit, 
                               y = fullFit, 
                               name = name)) + 
                    geom_point() + 
                    ggtitle(label = 'Comparing redundancy of Partial Models', subtitle = 'Percentage of Variance Explained by summation of partial and full')



additive_vs_full


sq %>% mutate(plus = conditionFit + strainFit) %>% ggplot(aes(x=plus, y=fullFit)) + geom_point() + geom_abline(intercept = 0, slope = 1, color = 'red')



parameters %>% filter(name == 'HSP12') %>% View()

parameters %>% dplyr::filter(abs(estimate)>1 & p.value < 0.01) %>% group_by(term, Model_Type) %>% summarise(num_genes = n())

parameters %>% dplyr::filter(abs(estimate)>1 & p.value < 0.01) %>% spread(Gene,term,estimate)

good_parameters = parameters %>% dplyr::filter(abs(estimate)>1 & p.value < 0.01)
lists = sapply(terms, function(x) select(filter(good_parameters,term==x),Gene))

estimates = parameters %>% select(Gene, term, estimate) %>% spread(term, estimate) %>% ungroup()
rownames(estimates) = estimates$Gene

estimates_matrix = as.matrix(select(estimates, -Gene))
col_mean = apply(estimates_matrix,2,mean)
#m = apply(parameter_matrix,1,function(x) x - col_m)
#m = t(m)
#heatmap3(m[m[,10]>1,-10], labRow = NA, scale='none')
heatmap3(parameter_matrix, labRow = NA, scale='none')


library(plotly)
library(d3heatmap)

shrink_large = function(data, max, min = -max){
data[data>max] = max
data[data<min] = min
return(data)
}

d3heatmap(a[str_length(rownames(a)) < 7,grep('Cond', colnames(a))], Rowv = NA, Colv=NA)

heatmaply(shrink_large(a[str_length(rownames(a)) < 7,
          grep('Cond', colnames(a))], 2),
          Rowv = NA,
          Colv=NA,
          scale_fill_gradient_fun= scale_fill_gradient2())
