library(sleuth)
library(dplyr)
library(gplots)
library(viridis)
library(ggplot2)
library(ggrepel)
library(UpSetR)

library("biomaRt")

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
mart
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)


t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Get all data
expression = sleuth_to_matrix(so, "obs_norm", "tpm")$data

# Get rid of genes that are lowly expressed
expression = expression[apply(expression, 1, min) > 0, ]

# Seperate out on Condition
hs_experiments = meta  %>%  filter(Condition == "HS")  %>%  select(sample)
wt_experiments = meta  %>%  filter(Condition == "NO", Drug == "NO")  %>%  select(sample)

expression_hs = expression[, unlist(hs_experiments)]
expression_wt = expression[, unlist(wt_experiments)]

# Remove bad sequencing run
expression_wt = expression_wt[, -which(colnames(expression_wt) %in% "RIM15_NO_NO")]

# Calculate foldchange
hs_fc = log2(expression_hs / expression_hs[, "WT_1NM_HS"])
wt_fc = log2(expression_wt / expression_wt[, "WT_NO_NO"])
#wt_fc = log2(expression_wt / apply(expression_wt, 1, mean)) # Could use this for mean
colnames(hs_fc) = unlist(lapply(colnames(hs_fc), function(x) strsplit(x, split="_")[[1]][1]))
colnames(wt_fc) = unlist(lapply(colnames(wt_fc), function(x) strsplit(x, split="_")[[1]][1]))

# Make heatmaps, remove wildtype becuase that was used for normalization
heatmap.2(cor(wt_fc[, -which(colnames(wt_fc) %in% "WT")]), col=viridis(100))
heatmap.2(cor(hs_fc[, -which(colnames(hs_fc) %in% "WT")]), col=viridis(100))





wtm = melt(wt_fc)
hsm = melt(hs_fc)
hsm$c = "HS"
wtm$c = "None"
all = rbind(hsm, wtm)
colnames(all) = c("Gene", "Strain", "Level", "Condition")
all = left_join(all, select(t2g, target_id, ext_gene), by=c("Gene"="target_id"))
all = all  %>%  dplyr::rename(Gene_Name=ext_gene)
ggplot(all) + geom_density(aes(x=Level,  fill = Condition))

hs = all  %>%  filter(Condition == "HS")
wt = all  %>%  filter(Condition == 'None')

# without heatshock
sampled_genes = rownames(expression)[sample(1:dim(expression)[1], 1, replace=F)]
  normal = ggplot() +
geom_dotplot(data = wt[wt$Gene %in% sampled_genes, ],  aes(x=Level),  alpha=0.3,  fill="blue") +
geom_density(data = wt[wt$Gene %in% sampled_genes, ],  aes(x=Level),  alpha=0.1,  colour="blue",  fill="blue") +
geom_text(data = hs[hs$Gene %in% sampled_genes, ],  aes(x=Level,  y=-.2,  label=Strain),  angle=90,  size=1.5,  check_overlap = TRUE) +
geom_point(data = hs[hs$Gene %in% sampled_genes, ],  aes(x=Level,  y=-.1)) +
coord_cartesian(ylim=c(-.3, 4)) +
facet_wrap(~Gene + Gene_Name)
normal
# with heatshock distribution
normal + geom_density(data = hs[hs$Gene %in% sampled_genes, ],  aes(x=Level),  alpha=0.3,  fill="red")

pdf("all.pdf")
for (i in 1:dim(expression)[1]) {
#for (i in 1:20) {
  sampled_genes = rownames(expression)[i]
  p = ggplot() +
  geom_dotplot(data = wt[wt$Gene %in% sampled_genes, ],  aes(x=Level),  alpha=0.3,  fill="blue") +
  geom_density(data = wt[wt$Gene %in% sampled_genes, ],  aes(x=Level),  alpha=0.1,  colour="blue",  fill="blue") +
  geom_text(data = hs[hs$Gene %in% sampled_genes, ],  aes(x=Level,  y=-.2,  label=Strain),  angle=90,  size=1.5,  check_overlap = TRUE) +
  geom_point(data = hs[hs$Gene %in% sampled_genes, ],  aes(x=Level,  y=-.1)) +
  coord_cartesian(ylim=c(-.3, 4)) +
  facet_wrap(~Gene + Gene_Name)
  print(p)
}
dev.off()



sampled_genes = rownames(expression)[sample(1:dim(expression)[1], 1, replace=F)]
ggplot(all[all$Gene %in% sampled_genes, ]) +
geom_density(aes(x = Level,  fill = Condition),  alpha = 0.3) +
coord_cartesian(xlim = c(-5,  5)) +
facet_wrap(~Gene + Gene_Name) +
geom_text(aes(x = Level,  y = 2,  label = Strain),  angle = 90,  size = 3) +
coord_cartesian(xlim = c(-2.5,  2.5))

ggplot(all) + geom_point(aes(x = Level,  fill = Strain), alpha = 0.25) + facet_wrap(~Condition)
ggplot(all) +
geom_density(aes(x = Level,  fill = Condition), alpha = 0.25) +
facet_wrap(~Strain) +
coord_cartesian(xlim = c(-5,  5))



expression_hs_log = log(expression_hs)
expression_hs_log = expression_hs_log[apply(expression_hs_log,  1,  function(x) all(is.finite(x))), ]
expression_hs_log_diff = expression_hs_log[, "WT_1NM_HS"]-expression_hs_log

sig = apply(expression_hs_log_diff, 2, function(x) ifelse(abs(x) > 2*sd(x), 1, -1))
# get rid of WT,  and TPK
sig = sig[, c(-1, -3)]
# Get rid of genes that are never significant
sig = sig[which(apply(sig, 1, sum) != -11), ]

km <- kmeans(sig,  50,  iter.max = 50)

sig = sig[order(km$cluster), ]


 image()



# melted_expression_hs = left_join(
#                          melt(expression_hs_log, value.name="expression"),
#                          melt(expression_hs_log_diff), by=c("Var1", "Var2"))  %>%
#                rename(transcript = Var1,
#                       condition = Var2,
#                       diff = value)  %>%
#                mutate(sig = diff > 1)

pdf("Heat Shock Expression.pdf")
for (i in 1:dim(expression_hs_log)[2]) {
  plot(expression_hs_log[, 1], expression_hs_log[, i],
    main = colnames(expression_hs_log)[i],
    xlab = "WT_1NM_HS",
    ylab = colnames(expression_hs_log)[i])
  idx = expression_hs_log_sig[, i]
  points(expression_hs_log[idx, 1], expression_hs_log[idx, i], col="red")
}
dev.off()

mmm = expression_hs_log_sig + 0
colnames(mmm) = unlist(lapply(colnames(mmm), function(x) strsplit(x, split="_")[[1]][1]))
upset(data.frame(mmm[, c(-1, -3)]), nsets=13, nintersects=40)

sig_genes = apply(mmm, 2, function(x) rownames(mmm[x == 1, ]))

lapply(seq_along(sig_genes),  function(ind,  list,  names) {
write(list[[ind]], paste(names[ind], "_ens_gene.txt", sep=""))
},  list = sig_genes,  names= names(sig_genes))

lapply(seq_along(sig_genes),  function(ind,  list,  names) {
write(apply(t2g[which(t2g[, 1]%in%sig_genes[[ind]]), ], 1, function(x) ifelse(nchar(x[3]) > 0, x[3], x[1])), paste(names[ind], "_ext_gene.txt", sep=""))
},  list = sig_genes,  names= names(sig_genes))
