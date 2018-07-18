library(tidyverse)
library(ggrepel)

mating_genes = read_lines('/Users/kieran/Thesis/kinase/input/genes_of_interest_and_gene_sets/GO_term_0019236_and_children_mating.tsv')
proteasome_genes = read_lines('/Users/kieran/Thesis/kinase/input/genes_of_interest_and_gene_sets/GO_term_0043161_and_children_proteosome.tsv')

output_path = '../../../../output/Images/figure_5'
dir.create(output_path)

# mating

genes %>% select(name, data,data_augment) %>%
  #filter(name %in% proteosome_module) %>%
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
#geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
coord_fixed(ratio=1) +
geom_abline(slope=1, intercept=0, alpha=.2) +
  geom_text_repel(data = . %>% filter(abs(.resid) > .2), color='black', ) + xlim(8.5,9.2) + ylim(8.5,9.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) -> p

ggsave(file.path(output_path, 'Figure_5a_scatter_mating.pdf'), plot = p)


# proteosome

genes %>% select(name, data,data_augment) %>%
  #filter(name %in% proteosome_module) %>%
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

ggsave(file.path(output_path, 'Figure_5a_scatter_proteosome.pdf'), plot = p)

#
# genes %>%
#   select(name, data, data_augment) %>%
#   #filter(name %in% proteosome_module) %>%
#   filter(name %in% mating_genes) %>%
#   unnest() %>%
#   #filter(Condition == 'Menadione') %>%
#   filter(Condition == 'Salt') %>%
#   mutate(Expression = Expression - mean(Expression[Strain == 'WT'])) %>%
#   ggplot(aes(x=Strain,
#              y = name,
#              fill = Expression)) +
#   geom_tile() +
#   #theme_Publication() +
#   xlab('Kianse') +
#   ylab('Gene') + scale_fill_gradient2() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
#
#   #geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
#   coord_fixed(ratio=1) +
#   geom_abline(slope=1, intercept=0, alpha=.2) +
#   geom_text_repel(data = . %>% filter(abs(.std.resid) > 3), color='black')
#
#
#
#
# genes %>% select(name, data,data_augment) %>%
#   filter(name %in% mating_genes) %>%
#   unnest() %>%
#   group_by(Condition, Strain_Code)%>%
#   summarize(Expression = mean(Expression), .fitted = mean(.fitted), .std.resid = mean(.std.resid)) %>%
#   ggplot(aes(x=Strain_Code,
#              y = Condition,
#              fill = .std.resid))  + geom_tile()
