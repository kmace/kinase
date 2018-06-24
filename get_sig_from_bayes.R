library(rethinking)
library(rstanarm)
library(shinystan)
library(tidyverse)
library(readr)

load('intermediate/images/normalized_data.RData')


sig_weight = function(name){
  f1 = vlog_blind[name,]
  df = data.frame(val = f1, Sample_Name = names(f1)) %>% left_join(meta)
  d = df %>% select(Condition, Strain, val)

  #simple = glimmer( val ~ Condition + Strain, data=d )
  interaction = glimmer( val ~ Condition + Strain + Condition:Strain, data=d )
  idx = (apply(interaction$d,2,var) != 0 & !grepl('CDC15', colnames(interaction$d)))
  fit <- stan_lm(val~., data = interaction$d[,idx], prior = R2(location = 0.7))

  #launch_shinystan(fit)

  #coeficients where every value is positive/negative
  #sign(apply(apply(as.data.frame(fit),2,range),2,prod)) %>% sort()
  sig = sign(apply(apply(as.data.frame(fit),2,range),2,prod))==1
  return(sig)
}

genes = names(which(apply(normalized_counts,1,mean)>35))

tab = purrr::map(genes, sig_weight)

mat = do.call(rbind, tab)

out_df = as.data.frame(mat, row.names = genes)

write.csv(x = out_df, file = 'true_interactions.csv')
