# ---- _1 ----
# General Data
t2g = load_transcripts_to_genes()
meta = load_meta_data()
# model_formula = ~ Strain + Drug + Condition
model_formula = ~ Deactivated_Kinase + Condition
# Slueth Data
sleuth_object = load_slueth_object(model.matrix(model_formula, data=meta),
                                   meta,
                                   t2g)
kallisto_counts = sleuth_to_matrix(sleuth_object,"obs_norm","est_counts")$data
# Deseq2 Data
star_counts = load_count_matrix(meta)
deseq_object = load_deseq_object(model_formula, meta, star_counts)
deseq_counts = counts(deseq_object, normalized=TRUE)
colnames(deseq_counts) = colnames(star_counts)
