library(dplyr)
library(tibble)


possible_motifs = expand.grid(Bg = c(0, 1),
                              Bs = c(0, 1),
                              Cg = c(-1, 0, 1),
                              Cs = c(-1, 0, 1),
                              Ds = c(-1, 0),
                              #Kg_on_E = c(-1, 0, 1),
                              #Ks_on_E = c(-1, 0, 1),
                              stringsAsFactors = F) %>% as_tibble #%>% group_by_at(vars(-K_type))

get_delta_E = function(df){
  General_Term = get_general_term(df$Bg, df$Cg)#, df$Kg_on_E)
  Specific_Term = get_specific_term(df$Bs, df$Cs, df$Ds)#, DF$Ks_on_E)
  equation = paste0('(', General_Term, ' * Kg_on_E) + (', Specific_Term, '* Kg_on_E)')
  return(equation)
}

  General_Term = get_general_term(Bg, Cg),#, Kg_on_E)
  Specific_Term = get_specific_term(Bs, Cs, Ds),#, Ks_on_E)
  equation = paste0('(', General_Term, ' * Kg_on_E) + (', Specific_Term, '* Kg_on_E)')



get_general_term = function(Bg, Cg){
if(Cg == 0){
  return(0)
} else {
  return(Bg - Cg)
}
}

get_specific_term = function(Bs, Cs, Ds){
  if(Cs == 0 & Ds == 0){
    return(0)
  } else if(Ds == 0){
    return(Bs - Cs)
  } else {
    return(Bs - Ds)
  }
}

possible_motifs %>% mutate(equation = map())

regulation_to_state = function(reg){
  if(reg=='+'){
    return('on')
  } else if(reg=='-'){
    return('off')
  } else{
    return(NA)
  }
}

regulation_to_state = Vectorize(regulation_to_state)

flip_reg = function(reg){
  if(reg=='+'){
    return('-')
  } else if(reg=='-'){
    return('+')
  } else if(reg=='0'){
    return('0')
  }
}

flip_reg = Vectorize(flip_reg)

reg_to_num = function(reg){
  if(reg=='+'){
    return(1)
  } else if(reg=='-'){
    return(0)
  } else if(reg=='0'){
    return(0)
  }
}

reg_to_num = Vectorize(reg_to_num)
possible_motifs %>%
  # Specify who may act:
  mutate(
    C_acts_on_G = C_on_G != '0',
    C_acts_on_K = C_on_K != '0',
    K_acts_on_G = K_on_G != '0') %>%
  # Account for Kinase modifications and tracking
  mutate(
    # Condition Specific Steps
    K_sense = if_else(C_acts_on_K, regulation_to_state(C_on_K), K_base), # Kinase state after Condition
    #C_changed_K = K_base != K_sense, # This is what we want to detect, when this is true
    # Analog Specific steps
    K_final = 'off') %>% #,
  #D_changed_K = K_sense != K_final,
  # Final checks
  #K_changed = K_base != K_final,
  #K_flipped = C_changed_K | D_changed_K) %>%
  # Figure out if something changed, and who was responsible for it
  mutate(
    C_acts_on_G_through_K_without_D = (K_base != K_sense) & K_acts_on_G,
    D_acts_on_G_through_K_without_C = K_base == 'on' & K_acts_on_G,
    D_blocks_C_action_on_G_through_K = C_acts_on_G_through_K_without_D & (K_sense == 'on'),
    D_and_C_are_redundant = C_acts_on_G_through_K_without_D & (K_sense == 'off')
  ) %>%
  # Account for changes in Gene Expression in the final state
  mutate(
    K_transitions = paste(K_base, K_sense, K_final, sep = '->')
  ) %>%
  mutate(
    synopsis = paste(paste('C has direct effect:', C_acts_on_G, sep = ' '),
                     paste('C direct effect sign: ', C_on_G, sep = ' '),
                     paste('D has a condition independent effect (K_base is on and affects G):', D_acts_on_G_through_K_without_C, sep = ' '),
                     paste('D effect sign (through K_final being off):', flip_reg(K_on_G), sep = ' '),
                     paste('C has an indirect effect through K (no D):', C_acts_on_G_through_K_without_D, sep = ' '),
                     paste('C indirect effect through K (no D) sign:', paste(C_on_K, K_on_G, sep='/'), sep = ' '),
                     paste('D blocks C effect:', D_blocks_C_action_on_G_through_K, sep = ' '),
                     paste('D and C are redundant:', D_and_C_are_redundant, sep = ' '),
                     paste('Kinase Transitions:', K_transitions, sep = ' '),
                     paste('Linear Expectation (C, K, R):', paste0('(', C_on_G,
                                                                   ', ', flip_reg(K_on_G),
                                                                   ', 0)'), sep = ' '),
                     paste('Complex Expectation (C, K, R):', paste0('(', C_on_G,
                                                                    ', ', flip_reg(K_on_G),
                                                                    ',', paste(C_on_K,
                                                                               K_on_G,
                                                                               sep='/') ,')'), sep = ' '),
                     sep = '\n')
    # g = purrr::map(., function(x) graph_from_data_frame(data.frame(from = c('c','c','k'),
    #                                      to = c('k', 'g', 'g'),
    #                                      strength = c(reg_to_num(x$C_on_K),
    #                                                   reg_to_num(x$C_on_G),
    #                                                   reg_to_num(x$K_on_G))),
    #                           directed=T,
    #                           vertices = nodes))
  ) -> outcomes
writeLines(outcomes$synopsis[2])
#synopsis = paste('Condition causes gene to go in a', C_on_G, 'direction.\n',
#					'With no analog, Kinase is in', Kinase_sensing_state, 'state, which causes the gene to go in a ', K_on_G, 'direciton when on.\n',
#			     'Kinase was basally', K_base_state, 'and',
#				 if_else(Condition_modified_state, 'was', 'was not'), 'flipped by the condition, and',
#			     if_else(Analog_modified_state, 'was', 'was not'), 'flipped by NMPP1 (to off).', sep=' ')
# NMPP1 either can break the basal level, or break condition signalling!
#)
