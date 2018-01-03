salt_genes = wt_results %>% dplyr::filter(condition=='Salt') %>% pull(name)
conditional_results %>%
  dplyr::filter(condition == 'Salt', name != 'HIS3') %>%
  dplyr::filter(name %in% salt_genes) %>%group_by(name) %>% mutate(n=n()) %>% arrange(desc(n)) %>% dplyr::filter(n<10) %>%
  select(Kinase, name) %>%
  split(.$Kinase) %>%
  lapply(function(x) pull(x,name)) -> sets

library(SuperExactTest)
res = supertest(sets, n = length(salt_genes), degree = c(2:4))
tab = summary(res)$Table %>% as_tibble() %>% arrange(P.value)  %>% mutate(log2FE = log2(FE))
tab %>% arrange(P.value)  %>% mutate(log2FE = log2(FE)) %>% `[`(1,7) %>% pbcopy()

good = res$P.value < 0.01 & res$overlap.sizes > 10
res$P.value = res$P.value[good]
res$overlap.sizes = res$overlap.sizes[good]
plot(res,degree=2, sort.by='p-value')


statistic = conditional_results %>% dplyr::filter(condition == 'Salt', Kinase == 'PBS2')
module.set = modules %>% split(.$module) %>% lapply(function(x) pull(x, name))
idx = ids2indices(module.set, statistic$name, remove.empty=TRUE)
cameraPR(statistic$stat, idx) %>% as_tibble(rownames = 'module') -> out
out
#spread(key=Kinase, value = pvalue) %>%
#as.data.frame() %>% remove_rownames() %>%
#column_to_rownames('name') %>% is.na() %>% `!` -> mat
