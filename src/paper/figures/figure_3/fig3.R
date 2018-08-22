# Figure 3
library(tidyverse)
library(dendextend)
library(circlize)
library(magrittr)
library(ComplexHeatmap)

per_condition_subset_results = read_csv('../../../../intermediate/per_condition_subset_results.csv')
wt_results = read_csv('../../../../intermediate/wt_results.csv')
t2g = read_csv('../../../../intermediate/t2g.csv')

get_cond = function(cond) {
  per_condition_subset_results %>%
    filter(condition == cond) %>%
    select(name, log2FoldChange, Kinase) %>%
    spread(value = log2FoldChange, key = Kinase)
}

get_kianse = function(kinase) {
  per_condition_subset_results %>%
    filter(Kinase == kinase) %>%
    select(name, log2FoldChange, condition) %>%
    spread(value = log2FoldChange, key = condition)
}

make_mat = function(x) {
  x %>%
    remove_rownames() %>%
    arrange(name) %>%
    column_to_rownames('name') %>%
    as.matrix()
}
My_Heatmap = function(x,
                      cluster_rows = F,
                      #col = mcol,
                      #colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                      col = colorRamp2(-5:5, divergent_colors),
                      ...){
  Heatmap(x,
          show_row_names = F,
          show_row_dend = F,
          col = col,
          cluster_rows = cluster_rows,
          clustering_method_columns = 'ward.D2',
          ...)
}

per_condition_subset_results %>%
  group_by(name) %>%
  summarise(num_missing  = sum(is.na(log2FoldChange))) %>%
  arrange(desc(num_missing)) %>%
  filter(num_missing == 0) %>%
  pull(name) -> good_genes

salt_genes = wt_results %>% dplyr::filter(condition=='Salt', pvalue < .1) %>% pull(name)
tm_genes = wt_results %>% dplyr::filter(condition=='Tunicamycin', pvalue < .1) %>% pull(name)

per_condition_subset_results %<>%
  filter(name %in% good_genes)

per_condition_subset_results_safe = per_condition_subset_results

## Tunic

per_condition_subset_results = per_condition_subset_results_safe %>% filter(name %in% tm_genes)

x = get_cond('Tunicamycin') %>% make_mat()
y = get_kianse('IRE1') %>% make_mat()
z = get_kianse('CDC15') %>% make_mat()
wt = wt_results %>% filter(condition == 'Tunicamycin') %>% select(name, log2FoldChange) %>% make_mat()
wt = wt[rownames(x), ]

tm_split = cutree(hclust(dist(cbind(x,y,z,wt)), method = 'ward.D2'),k=7)
tm_split %>% {write.csv(data.frame(name = names(.), cluster = .) %>% left_join(t2g), 'tm_clusters.csv')}

#Set output location:
output_path = '../../../../output/Images/figure_3'
dir.create(output_path)

png(file.path(output_path, 'Figure_3a_tunicamycin_heatmap.png'), width = 11000, height = 4000, res = 600)
My_Heatmap(x,
           name = 'Tunicamycin',
           split = tm_split) +
  My_Heatmap(y,
          name = 'IRE1') +
  My_Heatmap(z,
          name = 'CDC15') +
  My_Heatmap(wt,
          name = 'Wildtype')
dev.off()


## Salt


per_condition_subset_results = per_condition_subset_results_safe %>% filter(name %in% salt_genes)
x = get_cond('Salt') %>% make_mat()
y = get_kianse('PBS2') %>% make_mat()
z = get_kianse('TPK123') %>% make_mat()
wt = wt_results %>% filter(condition == 'Salt') %>% select(name, log2FoldChange) %>% make_mat()
wt = wt[rownames(x), ]

salt_split = cutree(hclust(dist(cbind(x,y,z,wt)), method = 'ward.D2'),k=6)
salt_split %>% {write.csv(data.frame(name = names(.), cluster = .) %>% left_join(t2g), 'salt_clusters.csv')}


png('Salt.png', width = 11000, height = 4000, res = 600)
My_Heatmap(x,
        name = 'Salt',
        split = salt_split) +
  My_Heatmap(y,
          name = 'PBS2'
          ) +
  My_Heatmap(z,
          name = 'PKA'
          ) +
  My_Heatmap(wt,
          name = 'Wildtype'
          )
dev.off()

sig_genes = per_condition_subset_results_safe %>% group_by(name) %>% summarise(mp = min(pvalue, na.rm = T)) %>% filter(mp<0.1) %>% pull(name)
per_condition_subset_results = per_condition_subset_results_safe %>% filter(name %in% good_genes) %>% filter(name %in% sig_genes)

pbs2 = get_kianse('PBS2') %>% make_mat()
pka = get_kianse('TPK123') %>% make_mat()
ssn3 = get_kianse('SSN3') %>% make_mat()
ypk1 = get_kianse('YPK1') %>% make_mat()
hog1 = get_kianse('HOG1') %>% make_mat()
cdc15 = get_kianse('CDC15') %>% make_mat()
ire1 = get_kianse('IRE1') %>% make_mat()
kss1 = get_kianse('KSS1') %>% make_mat()

k_mat = cbind(pbs2,
              pka,
              ssn3,
              ypk1,
              hog1,
              cdc15,
              ire1,
              kss1)

k_dist = dist(k_mat)
k_hc = hclust(k_dist, method = 'ward.D2')

kinase_split = cutree(k_hc,k=12)

kinase_split %>% {write.csv(data.frame(name = names(.), cluster = .) %>% left_join(t2g), 'kinase_clusters.csv')}

png('Kinase.png', width = 11000, height = 4000, res = 600)
My_Heatmap(pka,
           column_title = 'PKA',
           split = kinase_split, #14,
           #cluster_rows = as.dendrogram(k_hc),
           show_heatmap_legend = FALSE
           ) +
  My_Heatmap(pbs2,
             column_title = 'PBS2',
             show_heatmap_legend = FALSE
  ) +
  My_Heatmap(ypk1,
             column_title = 'YPK1',
             show_heatmap_legend = FALSE
  ) +
  My_Heatmap(ssn3,
             column_title = 'SSN3',
             show_heatmap_legend = FALSE
  ) +
  My_Heatmap(hog1,
             column_title = 'HOG1',
             show_heatmap_legend = FALSE
  ) +
  My_Heatmap(cdc15,
             column_title = 'CDC15',
             show_heatmap_legend = FALSE
  ) +
  My_Heatmap(ire1,
             column_title = 'IRE1',
             show_heatmap_legend = FALSE
  ) +
  My_Heatmap(kss1,
             column_title = 'KSS1',
             show_heatmap_legend = FALSE
  )
dev.off()

ref = genes %>% select(name, data, data_augment) %>% unnest() %>% select(name, Strain_Code, Condition, .std.resid) %>% rename(Kinase = Strain_Code, condition= Condition)

conditions %>% map(function(x){
  ref %>%
    filter(condition == x) %>%
    #filter(name %in% (wt_results %>% filter(condition == x & padj<0.05 & abs(log2FoldChange)>.5) %>% pull(name))) %>%
    select(name, Kinase, .std.resid) %>%
    spread(key = Kinase, value = .std.resid) %>%
    select(-name) %>%
    as.matrix() %>%
    cov()
}) -> out

names(out) = conditions

pdf('~/Desktop/epistasis.pdf')
map(conditions, function(x) {
Heatmap(out[[x]], col = colorRamp2(seq(-2,2,length.out = length(coolwarm_hcl)), coolwarm_hcl), column_title = x, clustering_method_rows = 'ward.D2',clustering_method_columns = 'ward.D2' )
})
dev.off()



genes %>%
select(name, data) %>%
unnest() %>%
group_by(name, Condition) %>%
mutate(DE = (Expression - mean(Expression[Strain == 'WT']))) %>%
ungroup() %>%
group_by(name) %>%
mutate(DE = DE/sd(DE)) %>%
select(name, Sample_Name, DE) %>%
spread(key = Sample_Name, value=DE) %>%
remove_rownames() %>%
column_to_rownames('name') %>%
as.matrix() -> exp_matrix2

exp_matrix2 = exp_matrix2[,sample_order]

tm_genes = tm_genes[tm_genes %in% rownames(exp_matrix2)]
x = exp_matrix2[tm_genes,((sm %>% filter(Condition == 'Tunicamycin') %>% pull(rowname)))]; colnames(x) = data.frame(rowname = colnames(x)) %>% left_join(sm) %>% pull(Kinase)
y = exp_matrix2[tm_genes,((sm %>% filter(Kinase == 'IRE1') %>% pull(rowname)))]; colnames(y) = data.frame(rowname = colnames(y)) %>% left_join(sm) %>% pull(Condition)
z = exp_matrix2[tm_genes,((sm %>% filter(Kinase == 'CDC15') %>% pull(rowname)))]; colnames(z) = data.frame(rowname = colnames(z)) %>% left_join(sm) %>% pull(Condition)

tm_split = cutree(hclust(dist(cbind(x,y,z)), method = 'ward.D2'),k=7)
as.data.frame(tm_split) %>% rownames_to_column('Gene') %>% arrange(tm_split) %>% write_csv(file.path(output_path, 'Figure_3a_tunicamycin_genes.csv'))

png(file.path(output_path, 'Figure_3a_tunicamycin_heatmap.png'), width = 11000, height = 4000, res = 600)
My_Heatmap(x,
name = 'Tunicamycin',
split = tm_split) +
My_Heatmap(y,
name = 'IRE1') +
My_Heatmap(z,
name = 'CDC15')
dev.off()

salt_genes = salt_genes[salt_genes %in% rownames(exp_matrix2)]
x = exp_matrix2[salt_genes,((sm %>% filter(Condition == 'Salt') %>% pull(rowname)))]; colnames(x) = data.frame(rowname = colnames(x)) %>% left_join(sm) %>% pull(Kinase)
y = exp_matrix2[salt_genes,((sm %>% filter(Kinase == 'HOG1') %>% pull(rowname)))]; colnames(y) = data.frame(rowname = colnames(y)) %>% left_join(sm) %>% pull(Condition)
z = exp_matrix2[salt_genes,((sm %>% filter(Kinase == 'PBS2') %>% pull(rowname)))]; colnames(z) = data.frame(rowname = colnames(z)) %>% left_join(sm) %>% pull(Condition)

na_split = cutree(hclust(dist(cbind(x,y,z)), method = 'ward.D2'),k=7)
as.data.frame(na_split) %>% rownames_to_column('Gene') %>% arrange(na_split) %>% write_csv(file.path(output_path, 'Figure_3b_salt_genes.csv'))

png(file.path(output_path, 'Figure_3b_salt_heatmap.png'), width = 11000, height = 4000, res = 600)
My_Heatmap(x,
           name = 'Salt',
           split = na_split) +
  My_Heatmap(y,
             name = 'HOG1') +
  My_Heatmap(z,
             name = 'PBS2')
dev.off()






hmap2 = Heatmap(exp_matrix2,
               name = '',
               cluster_columns = F,
               top_annotation = column_annotation,
               heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
               use_raster = FALSE,
               #col = colorRamp2(quantile_breaks(exp_matrix, 11), divergent_colors),
               col = colorRamp2(-5:5, divergent_colors),
               #split = split,
               show_row_names=F,
               show_column_names=F)

pdf(file.path(output_path, 'Figure_3c_contrast_heatmap.pdf'), width = 16, height = 9)
draw(hmap2, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

pdf(file.path(output_path, 'Figure_3c_skinny_contrast_heatmap.pdf'), width = 7, height = 9)
draw(hmap2, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

genes %>%
  select(name, data) %>%
  unnest() %>%
  group_by(name, Strain) %>%
  mutate(DE = (Expression - mean(Expression))) %>%
  ungroup() %>%
  group_by(name) %>%
  mutate(DE = DE/sd(DE)) %>%
  select(name, Sample_Name, DE) %>%
  spread(key = Sample_Name, value=DE) %>%
  remove_rownames() %>%
  column_to_rownames('name') %>%
  as.matrix() -> exp_matrix3

exp_matrix3 = exp_matrix3[,sample_order]

hmap3 = Heatmap(exp_matrix3,
                name = '',
                cluster_columns = F,
                top_annotation = column_annotation,
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
                use_raster = FALSE,
                #col = colorRamp2(quantile_breaks(exp_matrix, 11), divergent_colors),
                col = colorRamp2(-5:5, divergent_colors),
                #split = split,
                show_row_names=F,
                show_column_names=F)



pdf(file.path(output_path, 'Figure_3d_contrast_heatmap.pdf'), width = 16, height = 9)
draw(hmap3, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

pdf(file.path(output_path, 'Figure_3d_skinny_contrast_heatmap.pdf'), width = 7, height = 9)
draw(hmap3, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()


genes %>%
  select(name, data) %>%
  unnest() %>%
  group_by(name, Strain) %>%
  mutate(DE = (Expression - mean(Expression))) %>%
  ungroup() %>%
  group_by(name, Condition) %>%
  mutate(DE = (DE - mean(DE[Strain == 'WT']))) %>%
  ungroup() %>%
  group_by(name) %>%
  mutate(DE = DE/sd(DE)) %>%
  select(name, Sample_Name, DE) %>%
  spread(key = Sample_Name, value=DE) %>%
  remove_rownames() %>%
  column_to_rownames('name') %>%
  as.matrix() -> exp_matrix4

exp_matrix4 = exp_matrix4[,sample_order]

hmap4 = Heatmap(exp_matrix4,
                name = '',
                cluster_columns = F,
                top_annotation = column_annotation,
                heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
                use_raster = FALSE,
                #col = colorRamp2(quantile_breaks(exp_matrix, 11), divergent_colors),
                col = colorRamp2(-5:5, divergent_colors),
                #split = split,
                show_row_names=F,
                show_column_names=F)



pdf(file.path(output_path, 'Figure_3e_contrast_heatmap.pdf'), width = 16, height = 9)
draw(hmap4, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

pdf(file.path(output_path, 'Figure_3e_skinny_contrast_heatmap.pdf'), width = 7, height = 9)
draw(hmap4, heatmap_legend_side = "bottom")#, annotation_legend_side='bottom')
dev.off()

#### Figure 7

a = exp_matrix2[,((sm %>% filter(Kinase == 'TPK123') %>% pull(rowname)))]; colnames(a) = data.frame(rowname = colnames(a)) %>% left_join(sm) %>% pull(Condition)
b = exp_matrix2[,((sm %>% filter(Kinase == 'PBS2') %>% pull(rowname)))]; colnames(b) = data.frame(rowname = colnames(b)) %>% left_join(sm) %>% pull(Condition)

my_split = cutree(hclust(dist(cbind(a,b)), method = 'ward.D2'),k=7)
as.data.frame(my_split) %>% rownames_to_column('Gene') %>% arrange(my_split) %>% write_csv(file.path('../../../../output/Images/figure_7', 'Figure_7_genes.csv'))

png(file.path('../../../../output/Images/figure_7', 'Figure_7_anticor_heatmap.png'), width = 11000, height = 4000, res = 600)
My_Heatmap(a,
           name = 'PKA',
           split = my_split) +
  My_Heatmap(b,
             name = 'PBS2')
dev.off()

# Figure S3B
a = exp_matrix2[,((sm %>% filter(Kinase == 'YPK1') %>% pull(rowname)))]; colnames(a) = data.frame(rowname = colnames(a)) %>% left_join(sm) %>% pull(Condition)
b = exp_matrix2[,((sm %>% filter(Kinase == 'SSN3') %>% pull(rowname)))]; colnames(b) = data.frame(rowname = colnames(b)) %>% left_join(sm) %>% pull(Condition)

my_split = cutree(hclust(dist(cbind(a,b)), method = 'ward.D2'),k=7)
#as.data.frame(my_split) %>% rownames_to_column('Gene') %>% arrange(my_split) %>% write_csv(file.path('../../../../output/Images/figure_7', 'Figure_7_genes.csv'))

png(file.path('~/Desktop/Figure S3B.png'), width = 11000, height = 4000, res = 600)
My_Heatmap(a,
           name = 'YPK1',
           split = my_split) +
  My_Heatmap(b,
             name = 'SSN3')
dev.off()
