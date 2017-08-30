fs = full_std_resid_matrix
fs[abs(fs) < 2] = 0
fs_svd = svd(fs)
attach(fs_svd)

get_up = function(x, cut=2) {
norm = (x - mean(x)) / sd(x)
idx = which(norm > cut)
return(idx)
}

get_up(gene_out)
ups = apply(u,1,get_up)
ups$FUS1

get_top = function(x, n=5, high = TRUE, order_by_size = FALSE){
top_vals = head(sort(x, decreasing = high), n)
top_vals = x[x %in% top_vals]
if(order_by_size) {
top_vals = top_vals[order(abs(top_vals), decreasing = TRUE)]
}
return(top_vals)
}

get_top(u['FUS1',],10, high=FALSE)
get_top(u['AGA1',],10, high=FALSE)
get_top(u['KAR4',],10, high=FALSE)
get_top(u['FIG2',],10, high=FALSE)

gene = get_top(u[,8],30, F)
cond = get_top(v[,8],8, F)

std_residuals %>%
  filter(residual_type == 'Full') %>%
  left_join(t2g) %>%
  filter(name %in% names(gene)) %>%
  separate(col = sample_id, into = c('Cond', 'Strain'),sep = '_', remove = F) %>%
  group_by(Strain, Cond) %>%
  summarise(res = mean(residual)) %>%
  ggplot(aes(x=Strain,
             y=Cond,
             fill = res)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


get_top(u['MSC1',], 10, F)
get_top(u['CTT1',], 10, F)
get_top(u['DDR2',], 10, F)
get_top(u['TFS1',], 10, F)
get_top(u['GAD1',], 10, F)
get_top(u['SOL4',], 10, F)
get_top(u['PGM2',], 10, F)
get_top(u['DCS2',], 10, F)



#esr = t2g$name[match(scan('input/genes_of_interest_and_gene_sets/ESR/activated_ESR.txt', what = character()), t2g$Gene)]
#lapply(esr[esr %in% rownames(u)], function(x) get_top(u[x,], 10, F)) %>% lapply(names) %>% unlist %>% table() %>% sort


gene = get_top(u[,4],30)
cond = get_top(v[,4],8)


std_residuals %>%
  filter(residual_type == 'Full') %>%
  left_join(t2g) %>%
  filter(name %in% names(gene)) %>%
  separate(col = sample_id, into = c('Cond', 'Strain'),sep = '_', remove = F) %>%
  group_by(Strain, Cond) %>%
  summarise(res = mean(residual)) %>%
  ggplot(aes(x=Strain,
             y=Cond,
             fill = res)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
