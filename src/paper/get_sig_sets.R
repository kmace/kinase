library(tidyverse)
library(magrittr)
library(stringr)
library(DESeq2)
load('../../intermediate/images/normalized_data.RData')
throw = read_csv('../../intermediate/Samples_to_thow_out.csv') %>% left_join(meta)
dds_wt = dds_wt[,!(dds_wt$Sample_Name %in% throw$Sample_Name)]
dds_wt = dds_wt[,!grepl('WT4', dds_wt$Strain_Code)]

dds = dds[,!(dds$Sample_Name %in% throw$Sample_Name)]
dds = dds[,!grepl('WT4', dds$Strain_Code)]



wt_results = resultsNames(dds_wt)[-1] %>%
  map_dfr(~ lfcShrink(dds_wt, coef=.x) %>% # Use this when betaPrior was FALSE
            as.data.frame %>%
            as_tibble(rownames = 'name') %>%
            mutate(term = .x)) %>%
  mutate(condition = word(term,2,sep = '_'))

lapply(levels(dds$Condition)[-1], function(cond){
  dds_cond = dds[,dds$Condition==cond]
  dds_cond$Strain <- droplevels(dds_cond$Strain)
  dds_cond@design = ~Strain
  dds_cond = DESeq(dds_cond, parallel = T)
  res = lapply(resultsNames(dds_cond)[-1], function(x) {
    r = lfcShrink(dds_cond,
                  coef = x, parallel = TRUE) %>%
    #r = results(dds_cond, name = x, independentFiltering = FALSE) %>%
      as.data.frame() %>%
      as_tibble(rownames = 'name') %>%
      arrange(pvalue) %>%
      mutate(term = x, condition = cond)
  })
  results_tbl = do.call('rbind', res)
  return(results_tbl)
}) -> big_res
per_condition_subset_results = do.call('rbind', big_res) %>%
  mutate(Kinase = str_extract(term, '([A-z]{3}[0-9]{1,4})'))

write_csv(per_condition_subset_results, path='../../intermediate/per_condition_subset_results.csv')
write_csv(wt_results, path='../../intermediate/wt_results.csv')

