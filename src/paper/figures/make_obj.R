genes %>%
    select(name, data) %>%
    unnest() %>%
    group_by(name) %>%
    mutate(DE = (Expression - mean(Expression[Strain == 'WT' & Condition == 'YPD']))/sd(Expression)) %>%
    select(name, Sample_Name, DE) %>%
    spread(key = Sample_Name, value=DE) %>%
    remove_rownames() %>%
    column_to_rownames('name') %>%
    as.matrix() -> exp_matrix

sample_meta = meta[colnames(exp_matrix), ]
sample_meta$Kinase = sample_meta$Strain
sample_meta = as.data.frame(sample_meta) %>%
    select(Condition, Kinase, Sample_Name) %>%
    arrange(Condition, Kinase, Sample_Name) %>%
    remove_rownames() %>%
    column_to_rownames('Sample_Name')

sample_order = sample_meta %>% rownames()

exp_matrix = exp_matrix[,sample_order]
