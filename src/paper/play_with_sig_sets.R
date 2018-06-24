library(tidyverse)

# DESeq Data
per_condition_subset_results = read_csv('../../intermediate/per_condition_subset_results.csv')
wt_results = read_csv('../../intermediate/wt_results.csv')

# Set analysis
salt_genes = wt_results %>% dplyr::filter(condition=='Salt') %>% pull(name)
per_condition_subset_results %>%
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

source('figures/custom_smooth_function.R')

# Smooth by Condition
cond = 'Tunicamycin'

per_condition_subset_results %>%
#per_condition_subset_results %>% filter(name %in% (modules %>% filter(module == 'HOT1') %>% pull(name))) %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(condition == cond & wtp<0.1) %>%
  ggplot(aes(x=wt_change, y=log2FoldChange)) +
  stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
  geom_point() + geom_smooth(method = 'lm') + facet_wrap(~Kinase)

#for(cond in unique(wt_results$condition)){
#    for(tf in unique(factors$module)){
#        pdf(paste0('test/', tf,'_', cond,'.pdf'))
#        p = per_condition_subset_results %>% filter(name %in% (factors %>% filter(module == tf) %>% pull(name))) %>%


for(cond in levels(meta$Condition)[-1]){
  pdf(paste0(cond,'.pdf'))
  p = per_condition_subset_results %>%
    left_join(wt_results %>%
                dplyr::select(condition, log2FoldChange, name, padj) %>%
                dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
    dplyr::filter(condition == cond & wtp<0.1) %>%
    mutate(side = sign(wt_change)) %>%
    ggplot(aes(x=wt_change, y=log2FoldChange, group=side)) +
    stat_poly_eq(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "\n")),
                 formula = formula, parse = TRUE, coef.digits = 2, label.y.npc = "bottom") +
    geom_point() + geom_smooth(method='lm') + facet_wrap(~Kinase)
  print(p)
  dev.off()
}


# by Strain
kinase = 'CDC15'

per_condition_subset_results %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(Kinase == kinase & wtp<0.1 & abs(wt_change)>.5)%>%
  ggplot(aes(x=wt_change, y=log2FoldChange)) +
  stat_smooth_func(geom="text",hjust=0,parse=TRUE) +
  geom_point() + geom_smooth() + facet_wrap(~condition)


# by Strain
kinase = 'HOG1'
for(kinase in levels(meta$Strain)[-1]){
pdf(paste0(kinase,'.pdf'))
p = per_condition_subset_results %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(Kinase == kinase & wtp<0.1)%>%
  mutate(side = sign(wt_change)) %>%
  ggplot(aes(x=wt_change, y=log2FoldChange, group=side)) +
  stat_poly_eq(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "\n")),
                formula = formula, parse = TRUE, coef.digits = 2, label.y.npc = "bottom") +
  geom_point() + geom_smooth(method='lm') + facet_wrap(~condition)
print(p)
  dev.off()
}


genes_de_once = wt_results %>% filter(padj<0.01 & abs(log2FoldChange)>.5) %>% pull(name) %>% unique()
sample_meta = colData(dds_wt)
rld = rlog(dds_wt)
rlog = assay(rld)
rlog = rlog[genes_de_once, ]
source('../paper/figures/colors.R')


library(ComplexHeatmap)

sample_ha =  HeatmapAnnotation(sample_meta %>% as.data.frame %>% select(Condition),
                               col = master_col,
                               show_annotation_name = TRUE)

Heatmap(rlog - rowMeans(rlog[,sample_meta$Condition == 'YPD']),
        top_annotation = sample_ha,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        use_raster = F,
        show_column_names = F)

x = selectArea()
rownames(rlog)[x$row_order] %>% pbcopy


test = coef(dds_wt)
test = test[genes_de_once,]
colnames(test) %<>% word(2,sep = '_')
test = test[,-1]
Heatmap(test,
        #km = 5,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        use_raster = F)

x = selectArea()
rownames(test)[x$row_order] %>% pbcopy


library(modelr)

p = per_condition_subset_results %>%
  #per_condition_subset_results %>% filter(name %in% (modules %>% filter(module == 'HOT1') %>% pull(name))) %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(wtp<0.1 & condition == 'Tunicamycin' & Kinase %in% c('FUS3', 'IRE1', 'TPK123')) %>%
  ggplot(aes(x = wt_change, y = log2FoldChange)) + geom_point(size=.51) + geom_smooth(method='lm') + facet_wrap(~Kinase) + theme_bw() + ylim(c(-5,5))

ggsave(p,width = 15,height = 5, filename = '~/Desktop/slopes.pdf')

library(modelr)
per_condition_subset_results %>%
  #per_condition_subset_results %>% filter(name %in% (modules %>% filter(module == 'HOT1') %>% pull(name))) %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(wtp<0.1) %>% group_by(condition, Kinase) %>% nest() %>%
  mutate(mod = map(data,~lm(log2FoldChange ~ wt_change, data = .)),
         coefs = map(mod, ~coef(.)),
         slope = map_dbl(coefs, 'wt_change'),
         r2 = map2_dbl(mod, data, ~rsquare(.x,.y)),
         var = map_dbl(data, function(x) x %>% pull(log2FoldChange) %>% var)) -> fits

# No SDC
# Order cond by more hits (most to the right)
# Cluster Kinases
# Normalize by IRE1 TM
# 7 by 4 wide
# order by non abs mean columns
# get rid of gray background crosses.
library(magrittr)
fits %<>% filter(condition != 'SDC') %>%
  mutate(slope_norm = if_else(abs(slope)>abs(slope[Kinase == 'IRE1' & condition == 'Tunicamycin']),
                                      sign(slope) * abs(slope[Kinase == 'IRE1' & condition == 'Tunicamycin']),
                                      slope),
         var_norm = if_else(var>0.75*max(var),
                            0.75*max(var),
                            var))

fits %>% select(condition, Kinase, slope_norm) %>% spread(key = condition, value = slope_norm) %>% remove_rownames() %>% column_to_rownames('Kinase') %>% as.matrix() -> m
m %>% dist %>% hclust -> hc
Kinase_order = rownames(m)[hc$order]
#Condition_order = fits %>% group_by(condition) %>% summarize(m = mean((slope_norm[slope_norm<0]))) %>% arrange(desc(m))

fits$condition = factor(fits$condition, levels=c("Menadione",
                                                 "GlucoseDepletion",
                                                 "Rapamycin",
                                                 "Azetidine",
                                                 "Fluconazole",
                                                 "Salt",
                                                 "Tunicamycin",
                                                 "Heatshock"))
fits$Kinase = factor(fits$Kinase, levels=Kinase_order)

base = ggplot(fits, aes(y=Kinase, x=condition)) + #theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="black"),
        axis.line = element_line(colour = "black")
        )




# Fig 3 heatmaps
# try using plot_grid from cowplots, eg:
# weights = plot_grid(C,K, labels = c('B','C'), ncol=1)
# square_1 = plot_grid(expression, weights, resid, e_hat, labels=c('A','','D','E'), ncol=2)

per_condition_subset_results %>%
  group_by(name) %>%
  summarise(num_missing  = sum(is.na(log2FoldChange))) %>%
  arrange(desc(num_missing)) %>%
  filter(num_missing == 0) %>%
  pull(name) -> good_genes

per_condition_subset_results %<>%
  filter(name %in% good_genes)

library(dendsort)
cond_hm = function(C, filter = FALSE, col_dend = NULL, row_dend = NULL){
  per_condition_subset_results %>%
    filter(condition == C & name %in% good_genes) %>%
    select(name, log2FoldChange, Kinase) %>%
    spread(value = log2FoldChange, key = Kinase) %>%
    remove_rownames %>% column_to_rownames('name') %>%
    as.matrix %>% na.omit() -> mat
  if(filter){
    genes = wt_results %>% filter(condition == C & padj<0.1) %>% pull(name)
    genes_in = rownames(mat)[rownames(mat) %in% genes]
    mat = mat[genes_in,]
  }
  mat = mat[order(rownames(mat)),]
  if(is.null(col_dend)){
    col_dend = dendsort(hclust(dist(t(mat))))
  }
  if(is.null(row_dend)){
    row_dend = dendsort(hclust(dist((mat))))
  }

  p = Heatmap(mat,
              show_row_names = F,
              name = C,
              cluster_columns = col_dend,
              column_dend_reorder = FALSE,
              cluster_rows = row_dend,
              row_dend_reorder = FALSE)
  return(p)
}

kinase_hm = function(K){#, filter = FALSE){
  per_condition_subset_results %>%
    filter(Kinase == K) %>%
    select(name, log2FoldChange, condition) %>%
    spread(value = log2FoldChange, key = condition) %>%
    remove_rownames %>% column_to_rownames('name') %>%
    as.matrix %>% na.omit() -> mat
  # How would one filter on Kianse?
  #if(filter){
  #  mat = mat[wt_results %>% filter(condition == C & padj<0.1) %>% pull(name),]
  #}
  col_dend = dendsort(hclust(dist(t(mat))))

  p = Heatmap(mat,
              show_row_names = F,
              name = K,
              cluster_columns = col_dend,
              column_dend_reorder = FALSE)
  return(p)
}

library(cowplot)
library(ComplexHeatmap)

p = plot_grid(grid.grabExpr(draw(cond_hm('Salt',T))),
              grid.grabExpr(draw(cond_hm('Tunicamycin',T))),
              grid.grabExpr(draw(cond_hm('Heatshock',T))),
              grid.grabExpr(draw(kinase_hm('PBS2'))),
              grid.grabExpr(draw(kinase_hm('IRE1'))),
              grid.grabExpr(draw(kinase_hm('TPK123'))),
              labels=c('A','B','C','D','E','F'), nrow = 2)

pdf('~/Desktop/c_and_k_hmaps.pdf', width = 30, height = 20)
print(p)
dev.off()

p = plot_grid(grid.grabExpr(draw(cond_hm('Tunicamycin',T))),
              grid.grabExpr(draw(cond_hm('Salt',T))),
              grid.grabExpr(draw(kinase_hm('IRE1') + kinase_hm('CDC15'))),
              grid.grabExpr(draw(kinase_hm('PBS2') + kinase_hm('TPK123'))), nrow = 2)


salt_mat = per_condition_subset_results %>%
  filter(condition == 'Salt' & name %in% good_genes) %>%
  select(name, log2FoldChange, Kinase) %>%
  spread(value = log2FoldChange, key = Kinase) %>%
  remove_rownames %>% arrange(name) %>% column_to_rownames('name') %>%
  as.matrix %>% na.omit()
salt_mat = salt_mat[order(rownames(salt_mat)),]

row_dend_salt = dendsort(hclust(dist((salt_mat))))
col_dend_salt = dendsort(hclust(dist(t(salt_mat))))
dir.create('heatmaps_fig3')


for (filtered in c(TRUE, FALSE)) {
  for (row_order in c('onData', 'onSalt')) {
    for (col_order in c('onData', 'onSalt')) {
      for (condition in unique(per_condition_subset_results$condition)) {
        if(row_order == 'onData'){
          ro = NULL
        } else {
          ro = row_dend_salt
        }
        if(col_order == 'onData'){
          co = NULL
        } else {
          co = col_dend_salt
        }

        filename = paste0('heatmaps_fig3/',
                          condition, '_',
                          'filtered_', filtered, '_',
                          'genes_', row_order, '_',
                          'conditions_', col_order, '.pdf')


        pdf(filename)
        try(
        draw(cond_hm(condition, filtered, co, ro))
        )
        dev.off()
      }
    }
  }
}



mating = read_csv("~/Downloads/GOTerm_GeneOrganism-2.csv")
mating_t2g = t2g %>% filter(Gene %in% (mating %>% pull(Gene.secondaryIdentifier)))


g = (genes %>% select(name, data, data_augment) %>% unnest() %>% filter(name %in% mating_t2g$name) %>% group_by(name) %>% mutate(.fitted = .fitted - mean(Expression), Expression = Expression - mean(Expression)) %>% ungroup() %>% ggplot(aes(x = Expression, y=.fitted, color = Condition, label=Strain, text = name)) + geom_point() + geom_abline(slope=1, intercept=0, alpha=.2) + condition_color_scale)
ggsave('~/Desktop/mating_scatter.pdf')
ggp = ggplotly(g)
htmlwidgets::saveWidget(as_widget(ggp), "~/Desktop/mating_scatter.html")


g = (genes %>% select(name, data, data_augment) %>% unnest() %>% filter(name %in% mating_t2g$name) %>% group_by(name) %>% mutate(.fitted = .fitted - mean(Expression), Expression = Expression - mean(Expression)) %>% ungroup() %>% group_by(Condition, Strain) %>% summarize_all(mean) %>% ggplot(aes(x = Expression, y=.fitted, color = Condition, label=Strain, text = Strain)) + geom_point() + geom_abline(slope=1, intercept=0, alpha=.2) + geom_text_repel(data = . %>% filter(abs(.std.resid) > .51), color='black') + condition_color_scale)
ggsave('~/Desktop/mating_scatter_average.pdf')
ggp = ggplotly(g)
htmlwidgets::saveWidget(as_widget(ggp), "~/Desktop/mating_scatter_average.html")

protein_catabolism = read_csv("~/Downloads/GOTerm_GeneOrganism_ubiq_dep_prot_cat_proc.csv")
prot_t2g = t2g %>% filter(Gene %in% (protein_catabolism %>% pull(Gene.secondaryIdentifier)))

g = (genes %>% select(name, data, data_augment) %>% unnest() %>% filter(name %in% prot_t2g$name) %>% group_by(name) %>% mutate(.fitted = .fitted - mean(Expression), Expression = Expression - mean(Expression)) %>% ungroup() %>% ggplot(aes(x = Expression, y=.fitted, color = Condition, label=Strain, text = name)) + geom_point() + geom_abline(slope=1, intercept=0, alpha=.2) + condition_color_scale)
ggsave('~/Desktop/prot_cat_scatter.pdf')
ggp = ggplotly(g)
htmlwidgets::saveWidget(as_widget(ggp), "~/Desktop/prot_cat_scatter.html")


g = (genes %>% select(name, data, data_augment) %>% unnest() %>% filter(name %in% prot_t2g$name) %>% group_by(name) %>% mutate(.fitted = .fitted - mean(Expression), Expression = Expression - mean(Expression)) %>% ungroup() %>% group_by(Condition, Strain) %>% summarize_all(mean) %>% ggplot(aes(x = Expression, y=.fitted, color = Condition, label=Strain, text = Strain)) + geom_point() + geom_abline(slope=1, intercept=0, alpha=.2) + geom_text_repel(data = . %>% filter(abs(.std.resid) > 1), color='black') + condition_color_scale)
ggsave('~/Desktop/prot_cat_scatter_average.pdf')
ggp = ggplotly(g)
htmlwidgets::saveWidget(as_widget(ggp), "~/Desktop/prot_cat_scatter_average.html")

load('../../intermediate/images/paper_data.RData')

interactions = std_resid_matrix
interactions[abs(interactions)<2.5] = 0
interactions = sign(interactions)
library(philentropy)

jaccard_genes_d = distance(interactions, method = 'jaccard')
rownames(jaccard_genes_d) = rownames(std_resid_matrix)
colnames(jaccard_genes_d) = rownames(std_resid_matrix)

jaccard_conditions_d = distance(t(interactions), method = 'jaccard')
colnames(jaccard_conditions_d) = colnames(interactions)
rownames(jaccard_conditions_d) = colnames(interactions)

plot(hclust(as.dist(jaccard_conditions_d)))

library(ComplexHeatmap)
#Heatmap(interactions, show_row_names = F, show_column_names = F,
library(dendsort)


col_dend = dendsort(hclust(as.dist(jaccard_conditions_d)))
row_dend = dendsort(hclust(as.dist(jaccard_genes_d)))

p = Heatmap(interactions,
            show_column_names = F,
            #cluster_columns = col_dend,
            #column_dend_reorder = FALSE,
            show_row_names = F)#,
            #cluster_rows = row_dend,
            #row_dend_reorder = FALSE)


# write out the residual gene sets.
library(tidyverse)
dir.create('residual_sets')
std_residuals %>%
  filter(abs(residual)>2.5) %>%
  mutate(sign = sign(residual),
         label = paste(sample_id, if_else(sign == 1, 'underestimated', 'overestimated'), sep = '_')) %>%
  left_join(t2g) %>%
  split(.$label) %>%
  walk2(names(.), ~write_csv(.x, paste0('residual_sets/', .y, '.csv'))) %>%
  walk2(names(.), ~write(.x$Gene, paste0('residual_sets/', .y, '_proper.txt'))) %>%
  walk2(names(.), ~write(.x$name, paste0('residual_sets/', .y, '_common.txt')))

# Use python to upload these to SGD



go_res = read_csv('residual_sets/enrichments.csv') %>% select(-X1) %>% rename(pvalue = `p-value`) %>% separate(gene_set, into = c('Strain','Condition', 'Model_Estimate'), extra = 'drop', remove = F) %>% filter(pvalue<0.01)


# QQ plot

g = genes %>%
  select(data_augment) %>%
  unnest() %>%
  ggplot(aes(sample=.std.resid)) +
  stat_qq() +
  #geom_point() +
  geom_abline(slope=1, intercept=0, alpha=.2) + theme_classic()

ggsave(g, '~/Desktop/gg.png', dpi = 600)



t2g = read_csv('../../intermediate/t2g.csv')
iESR = t2g$name[t2g$Gene %in% scan('../../input/genes_of_interest_and_gene_sets/ESR/activated_ESR.txt', what = character())]
rESR = t2g$name[t2g$Gene %in% scan('../../input/genes_of_interest_and_gene_sets/ESR/repressed_ESR.txt', what = character())]

p = per_condition_subset_results %>%
  mutate(Group = if_else(Kinase %in% c('TPK123', 'PBS2'), Kinase, 'Other'),
         module = if_else(name %in% iESR, 'iESR', if_else(name %in% rESR, 'rESR','none'))) %>%
  filter(module != 'none') %>%
  ggplot(aes(x=log2FoldChange, group=Group, color = Group, fill = Group)) + geom_density(alpha=.2) + facet_wrap(~module) + theme_classic()

ggsave(p,'~/Desktop/pbs2_pka_hist.pdf')


per_condition_subset_results %>% filter(Kinase == 'PBS2' & name %in% genes$name) %>% select(log2FoldChange, condition, name) %>% spread(key = condition, value = log2FoldChange) %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix -> pbs2_mat
morpheus(na.omit(pbs2_mat))

library(morpheus)
per_condition_subset_results %>%
  dplyr::filter(Kinase == 'CDC15' & name %in% genes$name) %>%
  select(log2FoldChange, condition, name) %>%
  spread(key = condition, value = log2FoldChange) %>%
  remove_rownames() %>%
  column_to_rownames('name') %>%
  as.matrix -> cdc15_mat
morpheus(na.omit(cdc15_mat))



per_condition_subset_results %>%
  filter(condition == 'Salt' & name %in% genes$name & name %in% target_genes_na) %>%
  select(log2FoldChange, Kinase, name) %>%
  spread(key = Kinase, value = log2FoldChange) %>%
  remove_rownames() %>%
  column_to_rownames('name') %>%
  as.matrix -> na_mat
morpheus(na.omit(na_mat))

for(cond in unique(wt_results$condition)){
target_genes_cond = wt_results %>% filter(condition == cond & padj<0.1) %>% pull(name)

per_condition_subset_results %>%
filter(condition == cond & name %in% genes$name & name %in% target_genes_cond) %>%
select(log2FoldChange, Kinase, name) %>%
spread(key = Kinase, value = log2FoldChange) %>%
remove_rownames() -> cond_df
write_csv(cond_df, paste0('~/Desktop/expression_tables/', cond, '.csv'))
}


tf_files = list.files('processed_residual_sets/tfs/')

process_tf_file = function(fname){
  details = str_split(fname,pattern = '_')[[1]][1:3]
  tfs = read_lines(paste0('processed_residual_sets/tfs/',fname))
  if(length(tfs)>0){
  return(data.frame(Kinase = details[1], Condition = details[2], residual_type = details[3], tf = tfs))
  }
}

interactions = map_dfr(tf_files, process_tf_file)

process_tf_text_file = function(fname){
  details = str_split(fname,pattern = '_')[[1]][1:3]
  full_path = paste0('processed_residual_sets/meme/',fname, '/ame.txt')
  read_delim(full_path,
             delim = ' ',
             skip = 11,
             col_names = FALSE) ->tf
    if(dim(tf)[1]>0){
    tf = tf %>%
    select(c(7,8,10,13:15,20:22)) %>%
    rename(TF = X7,
           motif = X8,
           num_sequences_in_module_tested = X10,
           left_tail_pval = X13,
           right_tail_pval =  X14,
           two_tailed_pval = X15,
           corrected_left_tail_pval = X20,
           corrected_right_tail_pval =  X21,
           corrected_two_tailed_pval = X22)
      tf$Kinase = details[1]
      tf$Condition = details[2]
      tf$residual_type = details[3]
      return(tf)
    }
}

tf_files = list.files('processed_residual_sets/meme/')
interactions = map_dfr(tf_files, process_tf_text_file)

mes = genes %>% select(data_augment) %>% unnest()

png('~/Desktop/Residual_Plot_minimal.png')
mes %>%
  filter(Condition != 'Menadione') %>%
  ggplot(aes(x=.fitted,y=.std.resid, color = Condition)) +
  geom_point(size=0.1) +
  condition_color_scale +
  geom_abline(intercept = 2.5, slope = 0, alpha=.5) +
  geom_abline(intercept = -2.5, slope = 0, alpha=.5) +
  theme_minimal() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_blank(), legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank())
  #xlab("Predicted Level") +
  #ylab("Standardized Residual")
dev.off()



