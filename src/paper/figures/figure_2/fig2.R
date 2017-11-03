load('../../intermediate/images/paper_data.RData')
library(tidyverse)

weights = genes %>% select(name, weights) %>% unnest() %>% select(name, term, estimate)

K = weights %>%
    filter(grepl('Strain',term)) %>%
    spread(key = term, value = estimate)

C = weights %>%
    filter(grepl('Condition',term)) %>%
    spread(key = term, value = estimate)

M = genes %>%
    select(name, data, data_augment) %>%
    unnest()

E = M %>%
    select(name,
    Sample_Name,
    Differential_Expression) %>%
    spread(key = Sample_Name,
    value = Differential_Expression)

R = M %>%
    select(name,
    Sample_Name,
    .resid) %>%
    spread(key = Sample_Name,
    value = .resid)

gene_name = pull(E, name)


library(pheatmap)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
  breaks[!duplicated(breaks)]
}

my_p = function(mat, ...) {
    mat_breaks <- quantile_breaks(mat, n = 11)
    coolwarm_hcl <- colorspace::diverge_hcl(11,
                                            h = c(250, 10),
                                            c = 100,
                                            l = c(37, 88),
                                            power = c(0.7, 1.7))
    pheatmap(mat,
             color = coolwarm_hcl,
             scale = 'none',
             breaks = mat_breaks,
             border_color = NA,
             show_rownames = FALSE,
             show_colnames = TRUE,
             cluster_rows = F,
             cluster_cols = F,
             ...)
}

row_idx = hclust(dist(as.matrix(E[,-1])))$order

C = C[row_idx,]
E = E[row_idx,]
K = K[row_idx,]
R = R[row_idx,]

col_idx = hclust(dist(as.matrix(t(E[,-1]))))$order

c_names_ordered = colnames(E)[col_idx+1]

E = E[,c_names_ordered]
R = R[,c_names_ordered]
C = C[,-1]
K = K[,-1]



pdf('E.pdf')
my_p(as.matrix(E))
dev.off()
pdf('R.pdf')
my_p(as.matrix(R))
dev.off()
pdf('K.pdf')
my_p(as.matrix(K))
dev.off()
pdf('C.pdf')
my_p(as.matrix(C))
dev.off()
