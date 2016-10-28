tg = reg %>%
     filter(expert==TRUE, TF=='HSF1') %>%
     arrange(desc(log_prob_bind)) %>%
     select(Target) %>%
     head(40)

t2g$ext_gene[t2g$target_id %in% tg[,1]]
