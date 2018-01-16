library(tidyverse)
library(ggrepel)

mating_genes = read_csv('~/Desktop/mating_pathway.csv')
mating_genes = mating_genes$name

resid_modules %>% filter(module == 8) %>% pull(name) -> proteosome_module


genes %>% select(name, data,data_augment) %>%
  #filter(name %in% proteosome_module) %>%
  filter(name %in% mating_genes) %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .std.resid = mean(.std.resid)) %>%
    ggplot(aes(x=Expression,
               y = .fitted,
               color = Condition,
               label = Strain)) +
  geom_point(size=2) +
  condition_color_scale +
theme_Publication() +
xlab('Actual Expression') +
ylab('Predicted Expression') +
#geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
coord_fixed(ratio=1) +
geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.std.resid) > 1.5), color='black')

genes %>%
  select(name, data, data_augment) %>%
  #filter(name %in% proteosome_module) %>%
  filter(name %in% mating_genes) %>%
  unnest() %>%
  #filter(Condition == 'Menadione') %>%
  filter(Condition == 'Salt') %>%
  mutate(Expression = Expression - mean(Expression[Strain == 'WT'])) %>%
  ggplot(aes(x=Strain,
             y = name,
             fill = Expression)) +
  geom_tile() +
  #theme_Publication() +
  xlab('Kianse') +
  ylab('Gene') + scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


  #geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.std.resid) > 3), color='black')




genes %>% select(name, data,data_augment) %>%
  filter(name %in% mating_genes) %>%
  unnest() %>%
  group_by(Condition, Strain_Code)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .std.resid = mean(.std.resid)) %>%
  ggplot(aes(x=Strain_Code,
             y = Condition,
             fill = .std.resid))  + geom_tile()
