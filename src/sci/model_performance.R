load('../../input/images/model_parameters.RData')

mod_q = genes %>%
  select(Gene, name, data, ends_with('model')) %>%
  mutate_at(.vars=vars(ends_with('model')), 
            .funs = funs(glance = map(.,glance))) %>%
  select(Gene, name, ends_with('glance')) %>%
  gather(model_type, model_data, -Gene, -name) %>% 
  mutate(model_type = map(model_type, function(x) strsplit(x,'_')[[1]][1])) %>% 
  unnest()

# Distribution of explination of variance over genes for each non-base model
mod_q %>%
  filter(model_type != 'Base') %>%
  ggplot(aes(x = adj.r.squared, color=model_type, fill = model_type)) +
  geom_density(alpha=.5) + 
  ggtitle(label = 'Comparing fits of Models', subtitle = 'Percentage of Variance Explained by each model')


# Square up the r2
r2_sq = mod_q %>% ungroup() %>% select(name, model_type, adj.r.squared) %>%
  filter(model_type != 'Base') %>%
  spread(model_type, adj.r.squared)


# strain_vs_condition for each gene
ggplot(r2_sq,
       aes(x = Condition,
           y = Strain,
           name = name)) +
  geom_point() +
  ggtitle(label = 'Comparing fits of Partial Models', subtitle = 'Adjusted Percentage of Variance Explained by each model')

# condition_vs_full
ggplot(r2_sq,
       aes(x = Condition,
           y = Full,
           name = name)) +
  geom_point() +
  ggtitle(label = 'Improvement Full model has over Partial Models', subtitle = 'Adjusted Percentage of Variance Explained by each model') +
  geom_abline(intercept = 0, slope = 1, color = 'red')

# Strain_vs_full
ggplot(r2_sq,
       aes(x = Strain,
           y = Full,
           name = name)) +
  geom_point() +
  ggtitle(label = 'Improvement Full model has over Partial Models', subtitle = 'Adjusted Percentage of Variance Explained by each model') +
  geom_abline(intercept = 0, slope = 1, color = 'red')

# additive_vs_full - this is problematic! figure this out @TODO
# Works differently for r2 and adj.r2

r2_sq %>% 
  mutate(plus = Condition + Strain) %>% 
  ggplot(aes(x=plus, y=Full)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = 'red')


# Factor importance

getri = function(mod){
r = calc.relimp(mod, type = 'last')@last
return(r)
}


out = lapply(genes$Full_model, getri)
out
do.call(rbind, out)
lasts = do.call(rbind, out)
#pdf('')
plot(Strain~Condition, lasts, main='Additional R^2 from adding factor last')
