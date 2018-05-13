load('../../../../intermediate/images/paper_data.RData')
source('../make_obj.R')
source('../colors.R')
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

genes %>%
  select(Gene,weights) %>%
  unnest() %>%
  filter(grepl('Strain',term)) -> kinase_terms

kinase_terms = left_join(kinase_terms,t2g)


kinase_terms = kinase_terms %>%
    #filter(p.value < 0.05) %>%
    mutate(Kinase = gsub('Strain','',term))

ggplot(kinase_terms,
    aes(Kinase, fill = Kinase)) +
    geom_bar() #+
    #theme_Publication() +
    kinase_color_scale

kinase_term_matrix = kinase_terms %>% select(Kinase, name, estimate) %>% spread(key=Kinase,value=estimate)
kinase_term_matrix = kinase_term_matrix %>% remove_rownames() %>% column_to_rownames('name') %>% as.matrix()















Ma = M %>%
    filter(name=='AGA1') %>%
    select(Strain, Strain_Code, Condition, Differential_Expression, .resid)

Ea = Ma %>%
    select(-.resid) %>%
    spread(key = Condition, value = Differential_Expression) %>%
    select(-Strain_Code)
Ra = Ma %>%
    select(-Differential_Expression) %>%
    spread(key = Condition, value = .resid) %>%
    select(-Strain_Code)

measurement_col_ha = HeatmapAnnotation(data.frame(Condition = colnames(Ea)[-1]), col = master_col)
measurement_row_ha = rowAnnotation(data.frame(Kinase = Ea$Strain), col = master_col)

kinase = Ea$Strain

Ea = as.matrix(Ea[,-1])
Ra = as.matrix(Ra[,-1])

rep.row<-function(x,n){
   matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
   matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
Ca = C['AGA1',]
Ca['YPD'] = 0
Ca = rep.row(Ca, dim(Ea)[1])

Ka = K['AGA1',]
Ka['WT'] = 0
Ka = Ka[as.character(kinase)]

ea_hm = Heatmap(Ea, top_annotation=measurement_col_ha)
ra_hm = Heatmap(Ra, top_annotation=measurement_col_ha)
ca_hm = Heatmap(Ca, top_annotation=condition_ha)
ka_hm = Heatmap(Ka)

ea_hm + measurement_row_ha + ka_hm + ca_hm + ra_hm










# DESEq
library(DESeq2)
library(dplyr)
library(corrplot)
library(tidyverse)

load('../../../../intermediate/images/old_normalized_data.RData')
dds_salt = dds[,dds$Condition == 'Tunicamycin_YPD_Cocktail']
dds_salt$Strain <- droplevels(dds_salt$Strain)
dds_salt@design = ~ Strain
dds_salt = DESeq2::DESeq(dds_salt)

all_results = lapply(resultsNames(dds_salt),
function(x){
results(dds_salt, name=x, tidy=TRUE) %>%
as_tibble() %>%
left_join(t2g, by=c('row' = 'target_id')) %>%
arrange(padj) %>%
select(name, baseMean, log2FoldChange, padj) %>%
filter(padj<0.01)
}
)
names(all_results) = resultsNames(dds_salt)

do.call('rbind',
    lapply(seq_along(all_results),
    function(i) data.frame(test = names(all_results)[i],
                           all_results[[i]]))) %>%
  select(name, test, log2FoldChange) %>%
  spread(key=test, value=log2FoldChange, fill=0) %>%
  as_tibble() %>%
  select(-Intercept) %>%
  remove_rownames() %>%
  column_to_rownames('name') %>%
  cor() %>%
  corrplot(order='FPC')
