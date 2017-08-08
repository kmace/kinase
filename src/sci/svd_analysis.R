setwd("~/Thesis/kinase/src/shiny_apps/Individual_Gene_Explorer")
load('shiny.RData')
library(reshape2)
test = all %>% select(Gene_Name, Strain, Condition, Expression) %>% group_by(Gene_Name, Strain, Condition) %>% summarize(Expression = mean(Expression)) %>% ungroup()
td = acast(test, Gene_Name ~ Strain ~ Condition, value.var = 'Expression')

replace_na_with_mean = function(a){ 
    a[is.na(a)] = mean(a, na.rm = T); 
    #print(dim(a)); 
    return(a)
    }


for (i in 1:6692){ td[i,,] = replace_na_with_mean(td[i,,])}

# for aga1
single_gene = td['YIL117C',,] 
single_gene = single_gene - mean(single_gene)
out = svd(single_gene)
plot(1:length(out$d), out$d)
plot(1:length(out$d), log(out$d))

u = out$u
rownames(u) = rownames(single_gene)
d = out$d
v = out$v
rownames(v) = colnames(single_gene)

num_singular_values = 3
d[-(1:num_singular_values)] = 0
reconstruction = (u %*% diag(d) %*% t(v))
heatmap.2(reconstruction, col = viridis(20))
heatmap.2(single_gene, col = viridis(20))



barplot(names(td[1,,1]), apply(l1,1,mean))
barplot(apply(l1,1,mean))
barplot(apply(l1,1,mean), names.arg = names(td[1,,1]))



matplot(u, type = 'l')

aga1 = svd(clean['YNR044W',,])
aga1.svd = svd(aga1)
aga1
matplot(u)
barplot(u)
matplot(u)
matplot(u, type = 'l')
dim(u)
matplot(apply(u,1,function(col){col*diag(d)}), type = 'l')
u
d
aga1$u
aga1$u %*% aga1$d
plot(aga1$u %*% aga1$d)
matplot(aga1$u %*% aga1$d)
barplot(matplot(aga1$u %*% aga1$d))
matplot(aga1$u %*% aga1$d)
dim(aga1$u)
matplot(aga1$u %*% aga1$d)
plot(aga1$u %*% aga1$d)
barplot(aga1$u %*% aga1$d)
plot(aga1$u %*% aga1$d)
plot(aga1$d %*% t(aga1$v))


plot(aga1$d %*% t(aga1$v))

