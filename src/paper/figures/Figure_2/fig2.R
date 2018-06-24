library(tidyverse)

source('../custom_smooth_function.R')
source('../colors.R')

# Objects Required:
# per_condition_subset_results
# wt_results
#
#
#



# DESeq Data
per_condition_subset_results = read_csv('../../../../intermediate/per_condition_subset_results.csv')
wt_results = read_csv('../../../../intermediate/wt_results.csv')



#Set output location:
output_path = '../../../../output/Images/figure_2'
dir.create(output_path)




# Smooth by Condition
for (cond in c('Tunicamycin', 'Salt')) {
  p = per_condition_subset_results %>%
  #per_condition_subset_results %>% filter(name %in% (modules %>% filter(module == 'HOT1') %>% pull(name))) %>%
  left_join(wt_results %>%
              dplyr::select(condition, log2FoldChange, name, padj) %>%
              dplyr::rename(wt_change = log2FoldChange, wtp = padj)) %>%
  dplyr::filter(condition == cond & wtp<0.1) %>%
  ggplot(aes(x=wt_change, y=log2FoldChange)) +
  geom_point(size=.4) + geom_smooth(method = 'lm') + facet_wrap(~Kinase) + theme_classic() +
  xlab(paste0('WT-', cond, ' / WT-YPD')) +
  ylab(paste0('AS-', cond, ' / WT-', cond))

  ggsave(file.path(output_path, paste0(cond, '.pdf')), p)
}

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


pdf(file.path(output_path, 'slopes.pdf'), width = 7, height = 4.5)

# Full color slope
base +
  geom_tile(aes(fill = slope_norm)) +
  scale_fill_gradientn(na.value = 'black', colors = divergent_colors )
dev.off()

library(ggrepel)
pdf(file.path(output_path, 'r2_vs_slope.pdf'))
for(c in unique(fits$condition)){
  p = fits %>% filter(condition == c) %>% arrange(slope) %>% ggplot(aes(x=r2, y = slope, label=Kinase)) + geom_label_repel() + ggtitle(c)
  print(p)
}
dev.off()
