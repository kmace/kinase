library(stringr)
library(DESeq2)
library(tidyverse)
library(magrittr)
load('../../intermediate/images/normalized_data.RData')

res = lapply(resultsNames(dds_wt)[-1], function(x) {
  r = results(dds_wt,
              name = x,
              tidy = TRUE) %>%
    as_tibble %>%
    dplyr::rename(name = row) %>%
    arrange(pvalue) %>%
    dplyr::filter(padj<0.05) %>%
    mutate(term = x)
})
wt_results = do.call('rbind', res)
wt_results %<>% mutate(condition = word(term,2,sep = '_'))

lapply(levels(dds$Condition)[-1], function(cond){
  dds_cond = dds[,dds$Condition==cond]
  dds_cond$Strain <- droplevels(dds_cond$Strain)
  dds_cond@design = ~Strain
  dds_cond = DESeq(dds_cond, parallel = F)
  res = lapply(resultsNames(dds_cond)[-1], function(x) {
    r = results(dds_cond,
                name = x,
                tidy = TRUE) %>%
      as_tibble %>%
      dplyr::rename(name = row) %>%
      arrange(pvalue) %>%
      dplyr::filter(padj<0.05) %>%
      mutate(term = x, condition = cond)
  })
  results_tbl = do.call('rbind', res)
  return(results_tbl)
}) -> big_res
conditional_results = do.call('rbind', big_res)
conditional_results %<>% mutate(Kinase = str_extract(term, '([A-z]{3}[0-9]{1,4})'))

write_csv(conditional_results, path='../../intermediate/conditional_results.csv')
write_csv(wt_results, path='../../intermediate/wt_results.csv')


