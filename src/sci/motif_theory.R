library(dplyr)
library(tibble)

regulation = c('+','0', '-')
kinase_state = c('on', 'off')
kinase_type = c('WT', 'AS')
#NMPP1 = TRUE

possible_motifs = expand.grid(C_on_G = regulation,
							  C_on_K = regulation,
							  K_on_G = regulation,
							  K_type = kinase_type,
							  K_base = kinase_state,
							  stringsAsFactors = F) %>% as_tibble #%>% group_by_at(vars(-K_type))

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

possible_motifs %>%
  # Specify who may act:
  mutate(
     C_acts_on_G = C_on_G != '0',
	   C_acts_on_K = C_on_K != '0',
     D_acts_on_K = K_type == 'AS',  # Analog Sensitive
     K_acts_on_G = K_on_G != '0') %>%
  # Account for Kinase modifications and tracking
	mutate(
	   # Condition Specific Steps
	   K_sense = if_else(C_acts_on_K, regulation_to_state(C_on_K), K_base), # Kinase state after Condition
	   C_changed_K = K_base != K_sense, # This is what we want to detect, when this is true
	   # Analog Specific steps
	   K_final = if_else(D_acts_on_K, 'off', K_sense),
	   D_changed_K = K_sense != K_final,
	   # Final checks
	   K_changed = K_base != K_final) %>%
  # Figure out if something changed, and who was responsible for it
  mutate(
    C_acts_on_G_through_K = K_acts_on_G & C_changed_K,
    D_blocks_C_action_on_G_through_K = C_acts_on_G_through_K & D_changed_K,
    D_redundant_because_C = C_acts_on_G_through_K & !D_changed_K,
    D_acts_on_G_generally = D_acts_on_K & K_base == 'on' & K_acts_on_G
  ) %>%
  # Account for changes in Gene Expression in the final state
  mutate(
    C_contribution_to_G = C_on_G,
    K_contribution_to_G = if_else(K_final == 'on', K_on_G, flip_reg(K_on_G)),
    K_transitions = paste(K_base, K_sense, K_final, sep = '->'),
    G_movment = paste(C_contribution_to_G, K_contribution_to_G, sep='/')
    ) %>%
  mutate(
    synopsis = paste('C has direct effect:', C_acts_on_G,
                     'C has effect through K:', C_acts_on_G_through_K,
                     'D has general effect:', D_acts_on_G_generally,
                     'D reverses C\'s effect:', D_blocks_C_action_on_G_through_K,
                     'D redundant because C:', D_redundant_because_C,
                     'Kinase Transitions:', K_transitions,
                     sep = '\n')
  ) %>% View()

	   #synopsis = paste('Condition causes gene to go in a', C_on_G, 'direction.\n',
	   #					'With no analog, Kinase is in', Kinase_sensing_state, 'state, which causes the gene to go in a ', K_on_G, 'direciton when on.\n',
		#			     'Kinase was basally', K_base_state, 'and',
		#				 if_else(Condition_modified_state, 'was', 'was not'), 'flipped by the condition, and',
		#			     if_else(Analog_modified_state, 'was', 'was not'), 'flipped by NMPP1 (to off).', sep=' ')
	   # NMPP1 either can break the basal level, or break condition signalling!
   #)
