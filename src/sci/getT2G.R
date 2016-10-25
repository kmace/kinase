load_transcripts_to_genes = function(ds='scerevisiae_gene_ensembl') {
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = ds)
    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                         "ensembl_gene_id",
                                         "external_gene_name"), mart = mart)

    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    t2g <- t2g %>% dplyr::mutate(name=ifelse(ext_gene=='',ens_gene,ext_gene))
    return(t2g)
}

t2g_scerevisiae = load_transcripts_to_genes()
t2g_mmusculus = load_transcripts_to_genes('mmusculus_gene_ensembl')
