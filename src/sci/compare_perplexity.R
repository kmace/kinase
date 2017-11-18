library(tsne)

do_tsne = function(p){

  mat_meta = tibble(sample = colnames(mat)) %>%
    separate(sample, c('Condition', 'Strain'), "_", remove = FALSE)

  pca = prcomp(t(mat))
  pca = data.frame(pca$x, sample = rownames(pca$x))
  meta_pc = full_join(mat_meta, pca)
  t = tsne(meta_pc %>% select(starts_with('PC')), perplexity = p)
  colnames(t) = c('TSNE1', 'TSNE2')
  meta_pc = cbind(meta_pc, t)
  return(meta_pc)
}


#sample_tsne = do_tsne(40)


ps = c(seq(2,10,2), seq(15,60,5))
tsnes2 = lapply(tsnes,function(x) x %>% select(starts_with('TSNE'), Condition, Strain))
all = map2(ps, tsnes2, function(x,y) data.frame(perplexity = x, y))
all = do.call(rbind, all)

library(ggrepel)
ggplot(all,
       aes(x=TSNE1,
           y=TSNE2,
           color = Condition,
           label = Strain)) +
  geom_point() +
  #geom_text() +
  facet_wrap(~perplexity, scales = 'free')
