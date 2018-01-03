load('../../intermediate/images/normalized_data.RData')

library(ComplexHeatmap)
library(dendsort)
library(circlize)

for(target_condition in levels(dds_wt$Condition)[-1]){
pdf(paste0(target_condition, '.pdf'))
#target_condition = str_replace(target_condition,' ','.')
res = results(dds_wt, name = paste0("Condition_", target_condition, "_vs_YPD"))
rst_sub = rlog(dds[,dds$Condition == target_condition])
rlog = assay(rst_sub)
colnames(rlog) = colData(rst_sub)$Strain_Code
sub_median = apply(rlog[,grepl('WT', colnames(rlog))],1,median)
sub_genes = res %>% as.data.frame() %>% rownames_to_column('name') %>% filter(padj < .01 & abs(log2FoldChange) > 1) %>% left_join(t2g) %>% pull(name)
col_dend = dendsort(hclust(dist(t((rlog-sub_median)[sub_genes,]))))
draw(
  Heatmap((rlog- sub_median)[sub_genes,],
          show_row_names = F,
          name = 'Kinase vs. WT',
          cluster_columns = col_dend,
          column_dend_reorder = FALSE#,
          #col = colorRamp2(c(min((rlog- sub_median)[sub_genes,]),
          #                   0,
          #                   max((rlog- sub_median)[sub_genes,])), c("blue", "#EEEEEE", "red"))
          ) +
  Heatmap(res[sub_genes,]$log2FoldChange,
          show_row_names = F,
          name = paste0(target_condition, " vs. YPD")#,
          #col = colorRamp2(c(min(res[sub_genes,]$log2FoldChange),
          #                   0,
          #                   max(res[sub_genes,]$log2FoldChange)), c("blue", "#EEEEEE", "red"))
          )
)
dev.off()
}

