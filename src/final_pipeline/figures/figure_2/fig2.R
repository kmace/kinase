library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dendsort)
library(ggrepel)
library(cowplot)
load('../../../../intermediate/images/paper_data.RData')
source('../make_obj.R')
source('../colors.R')

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

E_hat = E - R


row_idx = dendsort(hclust(dist(E)))$order

C = C[row_idx,]
E = E[row_idx,]
E_hat = E_hat[row_idx,]
K = K[row_idx,]
R = R[row_idx,]

col_idx = dendsort(hclust(dist(t(E))))$order
c_names_ordered = colnames(E)[col_idx]

E = E[,c_names_ordered]
E_hat = E_hat[,c_names_ordered]
R = R[,c_names_ordered]

sample_meta = sample_meta[c_names_ordered,]


colnames(C) = gsub(pattern='Condition', replacement='', colnames(C))
colnames(K) = gsub(pattern='Strain', replacement='', colnames(K))


sample_ha =  HeatmapAnnotation(sample_meta,
                               col = master_col,
                               show_annotation_name = FALSE,
                               annotation_legend_param = list(Kinase = list(ncol = 1,
                                                                            title_position = "topleft",
                                                                            by_row=FALSE))
                                                                        )
condition_ha = HeatmapAnnotation(data.frame(Condition = colnames(C)), col = master_col)
kinase_ha = HeatmapAnnotation(data.frame(Kinase = colnames(K)), col = master_col)


quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm=T)
  breaks[!duplicated(breaks)]
}

# Consistnat coloring on E
col = colorRamp2(
  quantile_breaks(E, 11)[-c(1, 11)],
  coolwarm_hcl
)

# Consistnat coloring on E
test = quantile_breaks(abs(s), 5)
col = colorRamp2(
  c(-rev(test), 0, test)[-c(1, 11)],
  coolwarm_hcl
)


#col = colorRamp2(c(-3, 0, 3), c("green", "white", "red"))

make_hm = function(mat, ...){
    hm = Heatmap(mat,
        col = col,
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
    #column_title = expression(Delta~E[ij]),
    row_title = expression(gene[g]),
    top_annotation = sample_ha)

Heatmap(E,
        col = col,
        cluster_rows=T,
        cluster_columns=F,
        show_row_names=F,
        show_column_names=F,
        use_raster = TRUE,
        raster_quality = 5,
        name = 'Expression',
        #column_title = expression(Delta~E[ij]),
        row_title = expression(gene[g]),
        top_annotation = sample_ha)


ehat_hm = make_hm(E_hat,
    width=10, #unit(4, "in"),
    name = 'Expression',
    #column_title = expression(Delta~E[ij]),
    top_annotation = sample_ha)

r_hm = make_hm(R,
    width = 10, #unit(4, "in"),
    name = 'Residual',
    column_title = expression(R[ij]),
    top_annotation = sample_ha)

c_hm = make_hm(C,
    width = 4, #unit(4, "in"),
    name = 'Condition',
    #column_title = expression(C[i]),
    top_annotation = condition_ha)

k_hm = make_hm(K,
    width = 6, #unit(4, "in"),
    name = 'Kinase',
    #column_title = expression(K[j]),
    top_annotation = kinase_ha)

# png('heatmap.png', res = 500)
# e_hm + c_hm + k_hm + r_hm
# dev.off()

# png(filename='hm.png', width=24, height=5, units='in', res=300)
# e_hm + c_hm + k_hm + r_hm
# dev.off()

png(filename='hm.png', width=24, height=8, units='in', res=300)
#e_hm + c_hm + k_hm + r_hm
#e_hm + ehat_hm + r_hm
e_hm + ehat_hm + c_hm + k_hm

dev.off()

# png(filename='hm.png')
# e_hm + c_hm + k_hm + r_hm
# dev.off()

#
# theme_Publication <- function(base_size=14, base_family='') {
#       library(grid)
#       library(ggthemes)
#       (theme_foundation(base_size=base_size, base_family=base_family)
#        + theme(plot.title = element_text(face = "bold",
#                                          size = rel(1.2), hjust = 0.5),
#                text = element_text(),
#                panel.background = element_rect(colour = NA),
#                plot.background = element_rect(colour = NA),
#                panel.border = element_rect(colour = NA),
#                axis.title = element_text(face = "bold",size = rel(1)),
#                axis.title.y = element_text(angle=90,vjust =2),
#                axis.title.x = element_text(vjust = -0.2),
#                axis.text = element_text(),
#                axis.line = element_line(colour="black"),
#                axis.ticks = element_line(),
#                panel.grid.major = element_line(colour="#f0f0f0"),
#                panel.grid.minor = element_blank(),
#                legend.key = element_rect(colour = NA),
#                legend.position = "bottom",
#                legend.direction = "horizontal",
#                legend.key.size= unit(0.2, "cm"),
#                legend.spacing = unit(0, "cm"),
#                legend.title = element_text(face="italic"),
#                plot.margin=unit(c(10,5,5,5),"mm"),
#                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#                strip.text = element_text(face="bold")
#           ))
#
#}


# n = c('AGA1','AGA2','FUS1')
# n = iESR
# n = modules %>% select(module, genes) %>% unnest() %>% filter(module == 29) %>% pull(name)
#
# goTerms = readr::read_delim('~/geo_nlp_to_expression/GOTerm_GeneOrganism.tsv', delim = '\t', col_names = FALSE)
# colnames(goTerms) = c('sgd_gene_id', 'Gene', 'name', 'term_parent', 'term_parent_id', 'term', 'term_id', 'curation', 'verified', 'ontology', 'num', 'paper')
#
#
# n = goTerms %>% filter(grepl('shock', term, ignore.case = TRUE)) %>% pull(name) %>% unique()
#
# n = std_resid_matrix[(std_resid_matrix[,'SCH9_Menadione']) > 2,] %>% rownames()
#
# yeastnet = read_delim('~/Downloads/YeastNet.v3.txt', delim = '\t', col_names = c('Gene_1', 'Gene_2', 'weight'))
# yeastnet %>% na.omit() %>% graph.data.frame(directed = FALSE) -> g
# adj = as_adjacency_matrix(g, attr='weight')
fus1 = genes %>% filter(name == 'FUS1')


rng = range(c(unnest(fus1, data) %>% pull(Differential_Expression),
              unnest(fus1, data_augment) %>% pull(.resid)))
K = fus1 %>%
        select(weights) %>%
        unnest() %>%
        filter(grepl('Strain',term)) %>%
        mutate(Strain = gsub(pattern='Strain',replacement='', term)) %>%
        ggplot(aes(x=Strain, y = estimate, fill=Strain)) +
            geom_col() +
            theme_Publication() +
            strain_fill_scale +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, color = kinase_colors[-1]),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
            theme(legend.position="none")
C = fus1 %>%
        select(weights) %>%
        unnest() %>%
        filter(grepl('Condition',term)) %>%
        mutate(Condition = gsub(pattern='Condition',replacement='', term)) %>%
        ggplot(aes(x=Condition, y = estimate, fill=Condition)) +
            geom_col() +
            coord_flip() +
            theme_Publication() +
            condition_fill_scale +
            theme(axis.text.y = element_text(color = condition_colors[-1]),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
            theme(legend.position="none")

expression = fus1 %>%
        select(data) %>%
        unnest() %>%
        ggplot(aes(x = Strain, y=Condition, fill=Differential_Expression)) +
            geom_tile() +
            scale_fill_gradient2(low = coolwarm_hcl[1],
                                 mid = coolwarm_hcl[5],
                                 high = coolwarm_hcl[9],
                                 limits=c(rng[1],rng[2])) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, color = kinase_colors[-1]),
                  axis.text.y = element_text(color = condition_colors[-1]),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank()) +
            theme(legend.position="none")

resid = fus1 %>%
        select(data, data_augment) %>%
        unnest() %>%
        ggplot(aes(x = Strain, y=Condition, fill=.resid)) +
            geom_tile() +
            scale_fill_gradient2(low = coolwarm_hcl[1],
                                 mid = coolwarm_hcl[5],
                                 high = coolwarm_hcl[9],
                                 limits=c(rng[1],rng[2])) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, color = kinase_colors[-1]),
                  axis.text.y = element_text(color = condition_colors[-1]),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank()) +
            theme(legend.position="none")


e_hat = fus1 %>%
        select(data, data_augment) %>%
        unnest() %>%
        mutate(fitted_DE = Differential_Expression - .resid) %>%
        ggplot(aes(x = Strain, y=Condition, fill=fitted_DE)) +
            geom_tile() +
            scale_fill_gradient2(low = coolwarm_hcl[1],
                                 mid = coolwarm_hcl[5],
                                 high = coolwarm_hcl[9],
                                 limits=c(rng[1],rng[2])) +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, color = kinase_colors[-1]),
                  axis.text.y = element_text(color = condition_colors[-1]),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank()) +
            theme(legend.position="none")

weights = plot_grid(C,K, labels = c('B','C'), ncol=1)
square_1 = plot_grid(expression, weights, resid, e_hat, labels=c('A','','D','E'), ncol=2)
ggsave(square_1, 'square_1.pdf')

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

ggsave(p, filename = 'r2.pdf', height = 7, width = 8)



inlay = fus1 %>%
    mutate(r2 = map_dbl(performance, "r.squared")) %>%
    select(data, data_augment, r2) %>%
    unnest() %>%
    group_by(Condition, Strain) %>%
    summarise_all(mean) %>% ungroup() %>%
    ggplot(aes(x=Expression,y=.fitted, color = Condition, label = Strain)) +
        geom_point(size=2) +
        condition_color_scale +
        geom_text_repel(data = . %>% filter(abs(.std.resid) > 2), color='black') +
        theme_Publication() +
        xlab('Actual Expression') +
        ylab('Predicted Expression') +
        #geom_text(aes(mean(Expression), max(.fitted), label=paste("Gene: ", n, " | R2: ", round(r2, 3))), color = 'black') +
        coord_fixed(ratio=1) +
        geom_abline(slope=1, intercept=0, alpha=.2)

ggsave(inlay, filename = 'inlay.pdf', height = 8, width = 8)

# square_2 = ggdraw() +
#   draw_plot(inlay, 0  , 0  , 1  , 1  ) +
#   draw_plot(p,     0.5, 0.2, 0.3, 0.3) +
#   draw_plot_label(c("F", "G"),
#                   c(0, 1),
#                   c(.5, 0.5), size = 15)
#
# ggsave(square_2, )
