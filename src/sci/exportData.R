pinc = sleuth_to_matrix(so,"obs_norm","tpm")$data

control = pinc[,unlist(meta %>% filter(Drug=='NO',Condition == 'NO') %>% select(sample))]
just_drug = pinc[,unlist(meta %>% filter(Drug=='1NM',Condition == 'NO') %>% select(sample))]
drug_and_heat = pinc[,unlist(meta %>% filter(Drug=='1NM',Condition == 'HS') %>% select(sample))]

writeMat('expression.mat',
         control = control,
         control_conditions = colnames(control),
         control_transcript_id = rownames(control),
         just_drug = just_drug,
         just_drug_conditions = colnames(just_drug),
         just_drug_transcript_id = rownames(just_drug),
         drug_and_heat = drug_and_heat,
         drug_and_heat_conditions = colnames(drug_and_heat),
         drug_and_heat_transcript_id = rownames(drug_and_heat),
         transcript2gene = t2g[,c(1,3)])
