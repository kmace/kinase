library(rvest)
library(tidyverse)
library(data.table)


load_binary_yeastract_regulation_table = function(file, sep = ';') {
  reg = fread(file,
              header=T,
              sep=sep)
  # regulation = fread('../data/yeastract/RegulationMatrix_Documented_2013927.csv',
  #												 row.names=1,
  #												 quote='',
  #												 header=T,
  #												 sep=';')
  colnames(reg) = c('module', 'name')
  return (reg)
}

load_probBinding_yetfasco_regulation_table = function(file, use_expert=TRUE) {
  load(file)
  library(rvest)
  library(stringr)
  library(tidyr)
  txt = paste('<table>', str_c(datasetNames, collapse=' '), '</table>')
  html = read_html(txt)
  exp_table = html_table(html)[[1]]
  regT = data.frame(t(dataMat))
  regT[,'expert'] = FALSE
  regT[expert,'expert'] = TRUE
  regT[,'TF'] = exp_table[,2]
  regT[,'exp_name'] = exp_table[,1]
  reg = gather(regT, Target, log_prob_bind, -expert, -TF, -exp_name)
  return(reg)
}

load_pval_rickYong_regulation_table = function(file = '~/Desktop/Datasets/106_pvalbygene_ypd_v9.0.csv') {
  regT = read.csv(file,header=T, stringsAsFactors=F)
  reg = gather(regT,TF,pval, -ens_name,-ext_name,-description) %>%
    rename(Target=ext_name)
  reg$Target[reg$Target=="#REF!"] = reg$ens_name[reg$Target=="#REF!"]
  return(reg)
}
# Computational discovery of gene modules and regulatory networks
# Gifford lab


scrape_modules = function(p){
  YPD = read_html(p)
  YPD_cluster_names = YPD %>% html_nodes('h3') %>% html_text()
  YPD_cluster_cor = YPD %>% html_nodes(xpath = '/html/body/text()') %>% html_text()
  YPD_cluster_cor = YPD_cluster_cor[grepl('corr',YPD_cluster_cor)] %>% parse_number()
  YPD_cluster_membership = YPD %>% html_nodes('br+ table') %>% html_table(fill=TRUE)
  YPD_cluster_membership = lapply(YPD_cluster_membership, function(x){if(dim(x)[2]!=3){x = as.data.frame(matrix(ncol = 3))}; setNames(as_tibble(x), c('Gene', 'name', 'desc'))})
  YPD_cluster_meta = YPD %>% html_nodes('table+ table') %>% map(function(x) tryCatch(html_table(x, fill=TRUE),error=function(cond){return(data.frame(go=numeric(),term=character(),pval = numeric(), term_size=numeric(), module_term_size = numeric()))}))
  return(tibble(module=YPD_cluster_names,
                cor = YPD_cluster_cor,
                membership = YPD_cluster_membership,#,
                meta = YPD_cluster_meta
                ))
}

# YPD Modules
gifford_ypd_modules = scrape_modules('https://images.nature.com/original/nature-assets/nbt/journal/v21/n11/extref/nbt890-S3.htm')


# Rapamycin Modules
gifford_rap_modules = scrape_modules('https://images.nature.com/original/nature-assets/nbt/journal/v21/n11/extref/nbt890-S8.htm')



# segal modules
segal_modules = read_csv('input/external_datasets/module_definitions/segal_modules.csv')

yeastract_path = 'input/external_datasets/module_definitions/yeastract/'
yeastract = lapply(dir(yeastract_path, pattern = 'TwoColumn', full.names = T), load_binary_yeastract_regulation_table)
names(yeastract) = dir(yeastract_path, pattern = 'TwoColumn')

reg = load_probBinding_yetfasco_regulation_table('input/external_datasets/module_definitions/yetfasco/20120129_allMotifData1.02.rdat')
yetfasco = reg %>% filter(expert & log_prob_bind > -0.1)  %>% group_by(TF,Target) %>% summarise(ev_size = n()) %>% ungroup() %>% filter(str_length(TF)>0) %>% rename(module = TF, Gene = Target)



# 1) biochemical pathway common name 	- name of the biochemical pathway, as stored in SGD
# (mandatory)
# 2) enzyme name (optional)		- name of a specific enzyme (may be single or multi subunit)
# 3) E.C number of reaction (optional)	- Enzyme Commission identifier of the reaction, e.g. EC:1.1.1.1
# 4) gene name (optional)			- Gene name for the enzyme catalyzing the reaction, if identified
# 5) reference (optional)
biochem_pathway = read_tsv('input/external_datasets/module_definitions/biochemical_pathways.tab', col_names = F)
colnames(biochem_pathway) = c('module', 'enzyme', 'ec_num', 'name', 'ref')

pathways = read_tsv('input/external_datasets/module_definitions/Pathway_Genes.tsv', col_names = F)
