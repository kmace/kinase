library(tidyverse)
library(ggrepel)

mating_genes = read_lines('/Users/kieran/Thesis/kinase/input/genes_of_interest_and_gene_sets/GO_term_0019236_and_children_mating.tsv')
proteasome_genes = read_lines('/Users/kieran/Thesis/kinase/input/genes_of_interest_and_gene_sets/GO_term_0043161_and_children_proteosome.tsv')
heat_genes = read_lines('/Users/kieran/Thesis/kinase/input/genes_of_interest_and_gene_sets/GO_term_0009408_and_children_heat.tsv')
#iESR = t2g$name[t2g$target_id %in% scan('../../../../input/genes_of_interest_and_gene_sets/ESR/activated_ESR.txt', what = character())]
output_path = '../../../../output/Images/figure_5'
dir.create(output_path)

# hsp12
genes %>% select(name, data,data_augment) %>%
  filter(name %in% 'HSP12') %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .resid = Expression - .fitted) %>%
  ggplot(aes(x=Expression,
             y = .fitted,
             label = paste(Strain,'_',Condition, sep = ''))) +
  geom_point(size=2) +
  condition_color_scale +
  theme_Publication() +
  xlab('Average Actual Expression') +
  ylab('Average Predicted Expression') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > 2.8), color='black', ) + xlim(5,20.7) + ylim(5,20.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5a_scatter_hsp12.pdf'), plot = p)

# iESR
genes %>% select(name, data,data_augment) %>%
  filter(name %in% c('HSP12', 'PGM2', 'CTT1')) %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .resid = Expression - .fitted) %>%
  ggplot(aes(x=Expression,
             y = .fitted,
             label = paste(Strain,'_',Condition, sep = ''))) +
  geom_point(size=2) +
  condition_color_scale +
  theme_Publication() +
  xlab('Average Actual Expression') +
  ylab('Average Predicted Expression') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > 2.2), color='black', ) + xlim(7.7,17.6) + ylim(7.7,17.6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5b_scatter_msn2.pdf'), plot = p)

# aga1
genes %>% select(name, data,data_augment) %>%
  filter(name %in% 'AGA1') %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .resid = Expression - .fitted) %>%
  ggplot(aes(x=Expression,
             y = .fitted,
             label = paste(Strain,'_',Condition, sep = ''))) +
  geom_point(size=2) +
  condition_color_scale +
  theme_Publication() +
  xlab('Average Actual Expression') +
  ylab('Average Predicted Expression') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > 2.2), color='black', ) + xlim(7,13.5) + ylim(7,13.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5c_scatter_aga1.pdf'), plot = p)

# mating
genes %>% select(name, data,data_augment) %>%
  filter(name %in% mating_genes) %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .resid = Expression - .fitted) %>%
  ggplot(aes(x=Expression,
             y = .fitted,
             label = paste(Strain,'_',Condition, sep = ''))) +
  geom_point(size=2) +
  condition_color_scale +
  theme_Publication() +
  xlab('Average Actual Expression') +
  ylab('Average Predicted Expression') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > .2), color='black', ) + xlim(8.5,9.2) + ylim(8.5,9.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5d_scatter_mating.pdf'), plot = p)

# pre7
genes %>% select(name, data,data_augment) %>%
  filter(name %in% 'PRE7') %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .resid = Expression - .fitted) %>%
  ggplot(aes(x=Expression,
             y = .fitted,
             label = paste(Strain,'_',Condition, sep = ''))) +
  geom_point(size=2) +
  condition_color_scale +
  theme_Publication() +
  xlab('Average Actual Expression') +
  ylab('Average Predicted Expression') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > .7), color='black') + xlim(9.5,12.6) + ylim(9.5,12.6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5e_scatter_pre7.pdf'), plot = p)

# proteosome
genes %>% select(name, data,data_augment) %>%
  filter(name %in% proteasome_genes) %>%
  unnest() %>%
  group_by(Condition, Strain)%>%
  summarize(Expression = mean(Expression), .fitted = mean(.fitted), .resid = Expression - .fitted) %>%
  ggplot(aes(x=Expression,
             y = .fitted,
             label = paste(Strain,'_',Condition, sep = ''))) +
  geom_point(size=2) +
  condition_color_scale +
  theme_Publication() +
  xlab('Average Actual Expression') +
  ylab('Average Predicted Expression') +
  #geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > .3), color='black', ) + xlim(9.15,10.3) + ylim(9.15,10.3) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5f_scatter_proteosome.pdf'), plot = p)

