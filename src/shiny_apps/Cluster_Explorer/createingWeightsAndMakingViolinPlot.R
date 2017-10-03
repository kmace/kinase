load('../../../input/images/paper_data.RData')
ls()
genes
library(tidyverse)
models
library(broom)
#models %>% pull(Full_model) %>% map(tidy) -> weights
ls
weight
weights
?weighted.residuals
models %>% pull(Full_model) %>% map(tidy) -> weights
weights
models %>% mutate(weights = Full_model %>% map(tidy))
weights = models %>% mutate(weights = Full_model %>% map(tidy))
weights %>% select(weights)
w1 = weights %>% select(Gene, weights) %>% unnest()
w1
modules
clusters
w1 = weights %>% select(Gene, weights) %>% left_join(clusters)
w1
w1 %>% filter(cluster == 3)
w1 %>% filter(cluster == 3) %>% unnest()
w1 %>% filter(cluster == 3) %>% unnest() %>% group_by(term)
w1 %>% filter(cluster == 3) %>% unnest() %>% ggplot(aes(x = term, y = estimate) + geom_violin()
)
w1 %>% filter(cluster == 3) %>% unnest() %>% ggplot(aes(x = term, y = estimate)) + geom_violin()
#w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = ggplot(aes(x = term, y = estimate)) + geom_violin()
library(stringr)
?stringr::str_extract(
test = w1 %>% unnest() %>% pull(term)
test
str_extract(test, "[SC\(]")
str_extract(test, "[SC]")
str_extract(test, "[SC][a-z]")
str_extract(test, "[SC][a-z]*")
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin()
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) %>% ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin()
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) %>% ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin() + facet_grid(~term_type)
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) %>% ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin() + facet_grid(term_type~)
?facet_grid
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) %>% ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, )
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) %>% ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1)
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*")) %>% ggplot(aes(x = term, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free')
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type)) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free')
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free')
w1 %>% filter(cluster == 3) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
w1 %>% filter(cluster == 36) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
w1
w1$cluster %>% max
clusters
ls()
clusters
clusters$cluster %>% table
ls()
this_cluster
clusters_inspire
colnames(clusters_inspire)
w1
w2 = w1 %>% select(-cluster) %>% left_join(clusters_inspire)
w2
w1 %>% filter(Cluster == 36) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
w2 %>% filter(Cluster == 36) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
w2 %>% filter(name == 'ERO1')
w2 %>% filter(name == 'ERO1') %>% pull(Cluster)
w2 %>% filter(Cluster == 64) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
w2
w2 %>% rname(module = Cluster)
w2 %>% rename(module = Cluster)
w2 %>% rename(module = Cluster) %>% select(Gene, module, weights)
m = 36
w2 = w2 %>% rename(module = Cluster) %>% select(Gene, modul;5De, weights)
w2 = w2 %>% rename(module = Cluster) %>% select(Gene, modul, weights)
w2 = w2 %>% rename(module = Cluster) %>% select(Gene, module, weights)
w2
w2 %>% filter(Cluster == 64) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
w2 %>% filter(module == m) %>% unnest() %>% mutate(term_type = str_extract(term, "[SC][a-z]*"), term_name = str_replace(term, term_type, '')) %>% ggplot(aes(x = term_name, y = estimate, fill = term_type)) + geom_violin() + facet_wrap(~term_type, ncol = 1, scales='free') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
get_moduleWeightDistribution = function(m, weights){
  # For reference, this is the structure of weights since I havent made it yet
  # A tibble: 5,308 x 3
#       Gene module               weights
#      <chr>  <int>                <list>
#  1 YAL012W     81 <data.frame [38 x 5]>
#  2 YAL067C     39 <data.frame [38 x 5]>
#  3 YAL063C      5 <data.frame [38 x 5]>
#  4 YAL062W     62 <data.frame [38 x 5]>
#  5 YAL061W     73 <data.frame [38 x 5]>
#  6 YAL060W     70 <data.frame [38 x 5]>
#  7 YAL059W     67 <data.frame [38 x 5]>
#  8 YAL058W     16 <data.frame [38 x 5]>
#  9 YAL056W      1 <data.frame [38 x 5]>
# 10 YAL055W     44 <data.frame [38 x 5]>
# ... with 5,298 more rows
  weights %>%
    filter(module == m) %>%
    unnest() %>%
    mutate(term_type = str_extract(term, "[SC][a-z]*"),
           term_name = str_replace(term, term_type, '')) %>%
    ggplot(aes(x = term_name,
               y = estimate,
               fill = term_type)) +
    geom_violin() +
    facet_wrap(~term_type, ncol = 1, scales='free') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}
get_moduleWeightDistribution(22, w2)
?savehistory
?savehistory
savehistory(file = "createingWeightsAndMakingViolinPlot.R")
