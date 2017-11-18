# Module 85 -> independent Mating pathway
library(corrplot)
std_residuals %>% separate(sample_id, into = c('Kinase', 'Condition'), sep='_') %>% unite(col = Gene_Condition, name, Condition, remove = FALSE) %>% spread(key=Kinase, value = residual) -> std_residual_kinase_matrix
std_residual_kinase_matrix = left_join(std_residual_kinase_matrix, m)

std_residual_kinase_matrix %>% filter(mod == 85) %>% select(-mod) %>% select_if(is.numeric) %>% select_if(function(x) !all(is.na(x))) %>% cor(use = 'pairwise') %>% corrplot(order = 'AOE', diag = FALSE)

mgenes = genes %>% left_join(m) %>%
  filter(mod==77)

scatter = mgenes %>%
  #mutate(r2 = map_dbl(performance, "r.squared")) %>%
  select(data, data_augment) %>%
  unnest() %>%
  group_by(Condition, Strain) %>%
  summarise(.fitted = mean(.fitted), Expression = mean(Expression), .std.resid = mean(.std.resid)) %>%
  ungroup() %>%
  filter(Condition!='YPD' & Condition != 'Menadione') %>%
  ggplot(aes(x=Expression,y=.fitted, color = Condition, label = Strain)) +
  geom_point(size=2) +
  condition_color_scale +
  geom_text_repel(data = . %>% filter(abs(.std.resid) > 1), color='black') +
  theme_Publication() +
  xlab('Actual Expression') +
  ylab('Predicted Expression') +
  #geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
  coord_fixed(ratio=1) +
  geom_abline(slope=1, intercept=0, alpha=.2)

ggsave(scatter, filename = 'scatter.pdf')

strain = mgenes %>% select(weights) %>% unnest() %>% filter(grepl('Strain',term)) %>% mutate(Strain = gsub(pattern='Strain',replacement='', term)) %>% ggplot(aes(x=Strain, y = estimate, color=Strain)) + geom_boxplot() + theme_Publication() + strain_color_scale + theme(axis.text.x = element_text(angle = 90, hjust = 1, color = kinase_colors[-1])) + theme(legend.position="none")
ggsave(strain, filename = 'strain.pdf')

condition = mgenes %>% select(weights) %>% unnest() %>% filter(grepl('Condition',term)) %>% mutate(Condition = gsub(pattern='Condition',replacement='', term)) %>% ggplot(aes(x=Condition, y = estimate, color=Condition)) + geom_boxplot() + theme_Publication() + condition_color_scale + theme(axis.text.x = element_text(angle = 90, hjust = 1, color = condition_colors[-1])) + theme(legend.position="none")
ggsave(condition, filename = 'condition.pdf')
