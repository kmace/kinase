library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(RGraphics)
library(ggrepel)


load('../../input/images/normalized_data.RData')
source('../utils/0.1-data_loading_functions.R')
# If you want to look at all the data, then comment this out:
vsd = vsd[, colData(vsd)$Drug == 'Cocktail']

t2g = load_transcripts_to_genes()

colnames(vsd) = paste(vsd$Condition, vsd$Strain_Code, sep = "_")
out = assay(vsd)
out = as.data.frame(out)
out$Gene_Name = rownames(out)
all = tidyr::gather(out, Condition, Expression, -Gene_Name)
all = all %>% separate(Condition, c('Stress', 'Media', 'Drug', 'Strain'), "_")
all = all %>% dplyr::mutate(Condition = ifelse(Stress=='None',Media,Stress))
all$Name = t2g[match(all$Gene_Name, t2g$target_id),'name']
all$Description = t2g[match(all$Gene_Name, t2g$target_id),'description']

#is_outlier <- function(x) {
#  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
#}
is_outlier <- function(v, coef=1.5){
  quantiles <- quantile(v,probs=c(0.25,0.75))
  IQR <- quantiles[2]-quantiles[1]
  res <- v < (quantiles[1]-coef*IQR)|v > (quantiles[2]+coef*IQR)
  return(res)
}


all = all %>% group_by(Name, Condition) %>%
                mutate(is_condition_outlier = is_outlier(Expression),
                       condition_outlier = ifelse(is_condition_outlier,
                                                  as.character(Strain),
                                                  NA)) %>%
                ungroup() %>%
              group_by(Name, Strain) %>%
                mutate(is_strain_outlier = is_outlier(Expression),
                       strain_outlier = ifelse(is_strain_outlier,
                                               as.character(Condition),
                                               NA)) %>%
                ungroup()

all_data = all
reporters = read.table('../../input/reference/pathway_reporters.csv', sep='\t', header =T)
reporters$Gene = as.character(reporters$Gene)

make_plot = function(gene, pathway, all, t2g) {
dat = all %>% dplyr::filter(Name == gene)
p_heatmap = ggplot(dat, aes(y=Condition, x = Strain)) +
            geom_tile(aes(fill=Expression)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_cond_box = ggplot(dat, aes(x = Condition, y = Expression)) +
             geom_boxplot() +
             geom_text_repel(aes(label = condition_outlier), size=2, na.rm = TRUE) +
             coord_flip() +
             theme(axis.text.x = element_blank())
             

p_strain_box = ggplot(dat, aes(x = Strain, y = Expression)) +
               geom_boxplot() +
               geom_text_repel(aes(label = strain_outlier), size=2, na.rm = TRUE) +
               theme(axis.text.x = element_text(angle = 90, hjust = 1))


text = paste0(gene, " (", pathway, "): ", t2g[t2g$name == gene, "description"])
p_desc = ggplot() +
 annotate("text", x = 4, y = 25, size=8, label = text)

grid.arrange(p_heatmap + theme(legend.position="none", axis.text.x=element_blank(), axis.text.y=element_blank()),
          p_cond_box,
          p_strain_box + theme(legend.position="none", axis.text.y=element_blank()),
          textGrob(do.call(paste, c(as.list(strwrap(text, width = 0.7 * getOption("width"))), sep="\n"))), nrow=2, ncol=2)
}

output == FALSE
if(output){
pdf('../../output/reporter_genes.pdf', width=14, height=8)
for (i in 1:dim(reporters)[1]) {
    make_plot(reporters[i,2], reporters[i,1], all_data)
}
dev.off()
}