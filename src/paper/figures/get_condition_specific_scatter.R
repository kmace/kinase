iESR = t2g$name[t2g$Gene %in% scan('../../../../input/genes_of_interest_and_gene_sets/ESR/activated_ESR.txt', what = character())]
rESR = t2g$name[t2g$Gene %in% scan('../../../../input/genes_of_interest_and_gene_sets/ESR/repressed_ESR.txt', what = character())]


pdf('~/Desktop/relationship_esr.pdf', height = 9, width = 16)
values %>%
  group_by(name, Condition) %>%
  mutate(Differential_Expressoin = Expression - mean(Expression[Strain == 'WT'])) %>%
  filter(name %in% c(rESR, iESR)) %>%
  select(name, Strain_Code, Differential_Expression,Condition) %>%
  ungroup() %>%
  spread(key = Strain_Code, value = Differential_Expression) %>%
  mutate(Type = if_else(name %in% iESR, 'ESR Induced', 'ESR Repressed')) -> dat #%>%
  dat %>% ggplot(aes(x=HOG1,y=PBS2, color = Type)) +
  geom_point() +
  facet_wrap(~Condition, nrow = 2) +
  coord_fixed(ratio=1) +
  theme_Publication() +
  geom_abline(slope=1, intercept=0, alpha=.2)
dev.off()
