load('../../../intermediate/images/normalized_data.RData')

iESR = t2g$name[t2g$Gene %in% scan('../../../input/genes_of_interest_and_gene_sets/ESR/activated_ESR.txt', what = character())]
rESR = t2g$name[t2g$Gene %in% scan('../../../input/genes_of_interest_and_gene_sets/ESR/repressed_ESR.txt', what = character())]

normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column('name') %>%
  as_tibble() %>%
  gather(key = Sample_Name, value = norm_counts, -name) %>%
  mutate(Expression = log2(norm_counts + 1)) %>%
  left_join(meta) -> values


pdf('~/Desktop/relationship_esr.pdf', height = 9, width = 16)
values %>%
  group_by(name, Condition) %>%
  mutate(Differential_Expression = Expression - mean(Expression[Strain == 'WT'])) %>%
  ungroup() %>%
  select(name, Strain_Code, Differential_Expression, Condition) %>%
  dplyr::filter(name %in% c(rESR, iESR)) %>%
  spread(key = Strain_Code, value = Differential_Expression) %>%
  mutate(Type = if_else(name %in% iESR, 'ESR Induced', 'ESR Repressed')) -> dat #%>%
  dat %>% ggplot(aes(x=PBS2,y=TPK123, color = Type)) +
  geom_point() +
  facet_wrap(~Condition, nrow = 2) +
  coord_fixed(ratio=1) +
  theme_Publication() +
  geom_abline(slope=1, intercept=0, alpha=.2)
dev.off()
