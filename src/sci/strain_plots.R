vlog %>%
  as.data.frame() %>%
  rownames_to_column('target_id') %>%
  gather(key = Sample_Name, value = exp, -target_id) %>%
  left_join(meta) %>%
  group_by(Condition, target_id) %>%
  mutate(de = exp - median(exp[Strain == 'WT'])) %>%
  ungroup() -> diffE
for (target_strain in levels(dds$Strain)[-1]){
  diffE %>%
    filter(Strain == target_strain) %>%
    select(de, target_id, Condition) %>%
    spread(key = Condition, value = de) -> exp
  exp = as.matrix(exp[,-1])
  col_dend = dendsort(hclust(dist(t(exp))))
  pdf(paste0(target_strain, '.pdf'))
  draw(
    Heatmap(as.matrix(exp),
            show_row_names = F,
            name = 'Kinase vs. median(WT)',
            cluster_columns = col_dend,
            column_dend_reorder = FALSE,
            column_title = target_strain)
  )
  dev.off()
}

