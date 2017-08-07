make_tsne_plot = function(mat){
  
  mat_meta = tibble(cname = colnames(mat)) %>% 
    separate(cname, c('Condition', 'Strain'), "_", remove = FALSE)
  
  pca = prcomp(t(mat))
  pca = data.frame(pca$x, cname = rownames(pca$x))
  meta_pc = full_join(mat_meta, pca)
  t = tsne(meta_pc %>% select(starts_with('PC')))
  colnames(t) = c('TSNE1', 'TSNE2')
  meta_pc = cbind(meta_pc, t)
  ggplot(meta_pc, 
         aes(x=TSNE1, y=TSNE2, color=Condition, label = Strain)) + 
    geom_point(size=4) + 
    geom_text_repel(size = 1, 
                    color = 'black', 
                    point.padding = NA, 
                    box.padding = unit(0.01, "lines"))
  return(meta_pc)
}

pdf('tsnes.pdf')
base_dr = make_tsne_plot(base)
condition_dr = make_tsne_plot(condition)
strain_dr = make_tsne_plot(strain)
full_dr = make_tsne_plot(full)
dev.off()
