load('intermediate/images/normalized_data.RData')
library(DESeq2)
dds = dds[,colData(dds)$Strain != 'CDC15']
dds = dds[,colData(dds)$Condition != 'Menadione']
colData(dds)$Strain = droplevels(colData(dds)$Strain)
colData(dds)$Condition = droplevels(colData(dds)$Condition)
m1 = model.matrix(~Strain*Condition, colData(dds))
m2 = model.matrix(~Strain + Condition, colData(dds))

all.zero <- apply(m1, 2, function(x) all(x==0))
idx <- which(all.zero)
m1 <- m1[,-idx]
dds = DESeq(dds, test = 'LRT', full = m1, reduced = m2, betaPrior = FALSE)
res= results(dds)
save.image(file="ltr_test.RData")
res %>% as.data.frame() %>% rownames_to_column('target_id') %>% left_join(t2g) %>% filter(padj<0.05) %>% select(name, padj) -> res
res %>% filter(padj<0.05) %>% arrange(padj) %>% readr::write_csv(path='significant_lrt_genes.csv')
