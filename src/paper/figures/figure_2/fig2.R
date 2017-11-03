load('../../../../intermediate/images/paper_data.RData')
library(tidyverse)

weights = genes %>%
    select(name, weights) %>%
    unnest() %>%
    select(name, term, estimate)

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

C = C %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()
K = K %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()
R = R %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()
E = E %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()




row_idx = hclust(dist(E))$order

C = C[row_idx,]
E = E[row_idx,]
K = K[row_idx,]
R = R[row_idx,]

col_idx = hclust(dist(t(E)))$order
c_names_ordered = colnames(E)[col_idx]

E = E[,c_names_ordered]
R = R[,c_names_ordered]
sample_meta = as.data.frame(meta) %>% select(Strain, Condition) %>% dplyr::rename(Kinase = Strain)
sample_meta = sample_meta[c_names_ordered,]


colnames(C) = gsub(pattern='Condition', replacement='', colnames(C))
colnames(K) = gsub(pattern='Strain', replacement='', colnames(K))

library(RColorBrewer)
library(grDevices)

condition_color = RColorBrewer::brewer.pal(10,"Set3")
names(condition_color) = sort(unique(sample_meta$Condition))

kinase_color = colorRampPalette(brewer.pal(9, "Set1"))(29)
names(kinase_color) = sort(unique(sample_meta$Kinase))

master_col = list(Condition = condition_color,
                  Kinase = kinase_color)

sample_ha =  HeatmapAnnotation(sample_meta, col = master_col)
condition_ha = HeatmapAnnotation(data.frame(Condition = colnames(C)), , col = master_col)
kinase_ha = HeatmapAnnotation(data.frame(Kinase = colnames(K)), , col = master_col)


library(ComplexHeatmap)
library(circlize)
coolwarm_hcl <- colorspace::diverge_hcl(11,
                                        h = c(250, 10),
                                        c = 100,
                                        l = c(37, 88),
                                        power = c(0.7, 1.7))[-c(1, 11)]

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
  breaks[!duplicated(breaks)]
}

make_hm = function(mat, ...){
    hm = Heatmap(mat,
                 col = colorRamp2(
                            quantile_breaks(mat, 11)[-c(1, 11)],
                            coolwarm_hcl
                            ),
        cluster_rows=F,
        cluster_columns=F,
        show_row_names=F,
        show_column_names=F,
        use_raster = TRUE,
        raster_quality = 5,
        ...)
    return(hm)

}


e_hm = make_hm(E,
    width=10, #unit(4, "in"),
    name = 'Expression',
    column_title = expression(Delta~E[ij]),
    row_title = expression(gene[g]),
    top_annotation = sample_ha)

r_hm = make_hm(R,
    width = 10, #unit(4, "in"),
    name = 'Residual',
    column_title = expression(R[ij]),
    top_annotation = sample_ha)

c_hm = make_hm(C,
    width = 4, #unit(4, "in"),
    name = 'Condition',
    column_title = expression(C[i]),
    top_annotation = condition_ha)

k_hm = make_hm(K,
    width = 6, #unit(4, "in"),
    name = 'Kinase',
    column_title = expression(K[j]),
    top_annotation = kinase_ha)

# png('heatmap.png', res = 500)
# e_hm + c_hm + k_hm + r_hm
# dev.off()

# png(filename='hm.png', width=24, height=5, units='in', res=300)
# e_hm + c_hm + k_hm + r_hm
# dev.off()

png(filename='hm.png', width=24, height=8, units='in', res=300)
e_hm + c_hm + k_hm + r_hm
dev.off()

# png(filename='hm.png')
# e_hm + c_hm + k_hm + r_hm
# dev.off()


theme_Publication <- function(base_size=14, base_family='') {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))

}


p = genes %>%
    filter(Gene_mean > 6) %>%
    select(performance) %>%
    unnest() %>%
    ggplot(aes(x=r.squared)) +
    geom_density(fill='gray') +
    theme_Publication() +
    xlab(expression(paste(Individual~Model~Performance~"("~R^2~")"))) +
    ylab(expression('Density')) +
    scale_y_continuous(expand = c(0, 0))

ggsave(p, 'r2.pdf', height = 7, width = 8)

n = 'FUS1'



inlay = genes %>%
    filter(name==n) %>%
    mutate(r2 = map_dbl(performance, "r.squared")) %>%
    select(data, data_augment, r2) %>%
    unnest() %>%
    ggplot(aes(x=Expression,y=.fitted, color = Condition, label = Strain)) +
        geom_point(size=2) +
        condition_color_scale +
        geom_text_repel(data = . %>% filter(abs(.std.resid) > 2), color='black') +
        theme_Publication() +
        xlab('Actual Expression') +
        ylab('Predicted Expression') +
        geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
        coord_fixed(ratio=1) +
        geom_abline(slope=1, intercept=0, alpha=.2)

ggsave(inlay, 'inlay.pdf', height = 8, width = 8)
