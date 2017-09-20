library(dplyr)
library(tibble)

regulation = c('+','0', '-')
kinase_state = c('on', 'off')

possible_motifs = expand.grid(C_on_G = regulation,
							  C_on_K = regulation,
							  K_on_G = regulation,
							  K_base_state = kinase_state, stringsAsFactors = F) %>% as_tibble

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
mutate(# Condition Specific Steps
	   Condition_affects_Kinase = C_on_K != '0',
	   Kinase_sensing_state = if_else(Condition_affects_Kinase, regulation_to_state(C_on_K), K_base_state),
	   Condition_modified_state = K_base_state != Kinase_sensing_state,
	   # Analog Specific steps
	   Analog_affects_Kinase = TRUE,
	   Kinase_numb_state = 'off',
	   Analog_modified_state = Kinase_sensing_state != Kinase_numb_state,
	   # Possible outcomes
	   General_condition_effect = C_on_G,
	   Gener = 
	   synopsis = paste('Condition causes gene to go in a', C_on_G, 'direction.\n',
	   					'With no analog, Kinase is in', Kinase_sensing_state, 'state, which causes the gene to go in a ', K_on_G, 'direciton when on.\n',
					     'Kinase was basally', K_base_state, 'and',
						 if_else(Condition_modified_state, 'was', 'was not'), 'flipped by the condition, and',
					     if_else(Analog_modified_state, 'was', 'was not'), 'flipped by NMPP1 (to off).', sep=' ')
	   # NMPP1 either can break the basal level, or break condition signalling!
   )
