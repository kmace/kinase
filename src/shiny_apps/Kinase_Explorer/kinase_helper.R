library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

volcano_de = function(dds, contrast, t2g, pvalue_cutoff, fc_cutoff) {
  print(dds)
  print(contrast)
  print(pvalue_cutoff)
  print(fc_cutoff)
  res = results(dds, contrast=c("Strain",contrast,"WT"))
  res$Gene = t2g[match(rownames(res), t2g$target_id),]$name
  res$Description = t2g[match(rownames(res), t2g$target_id),]$description
  
  dat = as.data.frame(res) %>% 
    filter(!is.na(padj) & Gene != 'HIS3' & Gene != 'YOR203W') %>% 
    mutate(Significant = padj < pvalue_cutoff, 
           Large_Difference = abs(log2FoldChange) > fc_cutoff, 
           Sig_Diff = factor(paste(Significant, Large_Difference, sep="_"),
                             levels = c('FALSE_FALSE', 'FALSE_TRUE', 'TRUE_FALSE', 'TRUE_TRUE'),
                             labels = c('None', 'Large Difference', 'Significant', 'Large and Significant')))
  dat$Gene[dat$Sig_Diff != 'Large and Significant'] = NA
  
  if(table(dat$Sig_Diff)[4] > 200){
    dat$Gene = NA
  }
  
  plot = ggplot(dat,
                aes(
                  x = log2FoldChange,
                  y = -log10(padj),
                  color = Sig_Diff,
                  label = Gene
                )) + geom_point() +
    geom_text_repel(na.rm = TRUE, size=3) + 
    ggtitle(contrast)
  print(plot)
  return(plot)
}

volcano_tf = function(dds, contrast, t2g, score_metric, cutoff = NA, NA_rm = TRUE) {
  res = results(dds, contrast=c("Strain",contrast,"WT"))
  res$Gene = t2g[match(rownames(res), t2g$target_id),]$name
  res$Description = t2g[match(rownames(res), t2g$target_id),]$description
  
  dat = as.data.frame(res) %>% 
    filter(!is.na(padj) & Gene != 'HIS3' & Gene != 'YOR203W') 
  
  score_col = unlist(strsplit(score_metric,split='::'))[3]
  
  
  dat$Score = test[t2g[match(dat$Gene, t2g$name),]$target_id,score_col]
  
  if(!is.na(cutoff)){
    dat$Score = dat$Score > cutoff
  }
  
  if(!is.na(NA_rm)){
    dat = filter(dat, !is.na(Score))
  }
  
  plot = ggplot(dat,
                aes(
                  x = log2FoldChange,
                  y = -log10(padj),
                  color = Score,
                  label = Gene
                )) + geom_point() + 
    ggtitle(contrast)
  print(plot)
  return(plot)
}