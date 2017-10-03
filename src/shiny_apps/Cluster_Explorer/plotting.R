get_moduleSummary = function(m, modules, tfs){
  size = modules %>%
    filter(module == m) %>%
    nrow()
  tf = tfs %>%
    filter(module == m) %>%
    pull(TF)

  size_text = paste0('Module Size: ', size, ' Genes')
  tfs_text = paste0('TFs: ', paste(tf, collapse = "; "))

  summary_text = paste0(size_text,
                        '\n',
                        tfs_text)

  return(summary_text)
}

get_moduleSimilarity = function(m, MEsCor){
  mod_cor = MEsCor[, paste0('ME', as.character(m))]
  top_10 = rev(names(head(tail(sort(abs(mod_cor)),11),10)))
  top_cor = mod_cor[top_10]
  tibble(Module = names(top_cor),
         Corrolation = top_cor) %>%
    mutate(sign = if_else(sign(Corrolation) == 1, 'Pos', 'Neg')) %>%
    ggplot(aes(x = reorder(Module, abs(Corrolation)),
               y = Corrolation,
               fill=sign)) +
    geom_bar(stat='identity') +
    theme(legend.position="none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) +
    coord_flip()
}

get_moduleWeightDistribution = function(m, weights){
  # For reference, this is the structure of weights since I havent made it yet
  # A tibble: 5,308 x 3
#       Gene module               weights
#      <chr>  <int>                <list>
#  1 YAL012W     81 <data.frame [38 x 5]>
#  2 YAL067C     39 <data.frame [38 x 5]>
#  3 YAL063C      5 <data.frame [38 x 5]>
#  4 YAL062W     62 <data.frame [38 x 5]>
#  5 YAL061W     73 <data.frame [38 x 5]>
#  6 YAL060W     70 <data.frame [38 x 5]>
#  7 YAL059W     67 <data.frame [38 x 5]>
#  8 YAL058W     16 <data.frame [38 x 5]>
#  9 YAL056W      1 <data.frame [38 x 5]>
# 10 YAL055W     44 <data.frame [38 x 5]>
# ... with 5,298 more rows
  weights %>%
    filter(module == m) %>%
    unnest() %>%
    mutate(term_type = str_extract(term, "[SC][a-z]*"),
           term_name = str_replace(term, term_type, '')) %>%
    ggplot(aes(x = term_name,
               y = estimate,
               fill = term_type)) +
    geom_violin() +
    facet_wrap(~term_type, ncol = 1, scales='free') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

get_moduleExpression = function(m, module_raw_residuals){
  module_raw_residuals %>%
    filter(residual_type == 'Base', module == m) %>%
    separate(sample_id, c('Condition', 'Kinase'), '_') %>%
    ggplot(aes(x = Kinase,
               y = Condition,
               fill = avg_residual)) +
    geom_tile() +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = "log2") +
    ggtitle("Differnetial Gene Expression")
}

get_moduleResidual = function(m, module_std_residuals){
  module_std_residuals %>%
    filter(residual_type == 'Full', module == m) %>%
    separate(sample_id, c('Condition', 'Kinase'), '_') %>%
    ggplot(aes(x = Kinase,
               y = Condition,
               fill = avg_residual)) +
    geom_tile() +
    scale_fill_gradient2(low = 'cyan', mid = 'black', high = 'yellow') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(fill = expression(sigma)) +
    ggtitle("Residuals error in full model")
}

get_membershipTable = function(m, modules) {
  return(modules %>%
           filter(module == m) %>%
           select(name, Gene, description) %>%
           rename(Common_Name = name,
                  Description = description))
}

get_GOEnrichmentTable = function(m, tab){
  tab %>%
    filter(module == m) %>%
    rename(Num_Module_Genes_in_Term = nModGenesInTerm,
           Num_Genes_in_Term = bkgrTermSize,
           Bonferoni_pval = BonferoniP,
           GO_Term = termName) %>%
    select(GO_Term,
           Num_Module_Genes_in_Term,
           Num_Genes_in_Term,
           Bonferoni_pval)
}

get_memePage = function(m){
  path = paste0('WGCNA_modules/meme_output/cluster_', m, '_genes/ame.html')
  return(tags$iframe(src=path, height=600, width="100%"))
}
