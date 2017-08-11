load('../../input/images/model_parameters.RData')

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

# For the samples:
library(tsne) # original implementation, exact tsne, slow
pdf('tsnes.pdf')
base_dr = make_tsne_plot(base)
condition_dr = make_tsne_plot(condition)
strain_dr = make_tsne_plot(strain)
full_dr = make_tsne_plot(full)
dev.off()


#for the genes:
library(Rtsne) # much faster, but aproximate
# We will let Rtsne do the PCA for us.
cols = viridis::viridis(2)
out = Rtsne(base)
colnames(out$Y) = c('TSNE1', 'TSNE2')
gene_cols = cols[rownames(base) %in% 
                   scan('../../input/genes_of_interest_and_gene_sets/ESR/activated_ESR.txt', 
                        what = character()) + 1]
plot(TSNE1~TSNE2, 
     data = out$Y, 
     col = gene_cols, 
     main = 'activated ESR')

gene_cols = cols[rownames(base) %in% 
                   scan('../../input/genes_of_interest_and_gene_sets/ESR/repressed_ESR.txt', 
                        what = character()) + 1]
plot(TSNE1~TSNE2, 
     data = out$Y, 
     col = gene_cols,
     main = 'repressed ESR')

