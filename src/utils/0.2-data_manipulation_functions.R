get_gene_list_from_reg_table = function(reg_table, tf, cutoff=NULL) {
    tf_rows = which(reg_table$TF == tf)
    subset = reg_table[tf_rows,]
    if(!is.null(cutoff)){
        good_rows = which(subset$pval <= cutoff)
        subset = subset[good_rows, ]
    }
    return(subset$Target)
}

add_reg_metadata = function(reg, t2g, expression){
    #t2g <- t2g %>% dplyr::mutate(name=ifelse(ext_gene=='',ens_gene,ext_gene))
    reg$index = match(reg$Target,t2g$name)
    reg$name = t2g$target_id[reg$index]
    reg$expression_index = match(reg$name,rownames(expression))
    return(reg)
}


filter_undetected_genes = function(data) {
    no_zeros = apply(data,1,min) > 0
    return(data[no_zeros,])
}
# ---- _2 ----
filter_and_foldchange = function(data) {
    filtered = filter_undetected_genes(data)
    fold_change = apply(filtered, 1, function(x) log2(x/mean(x)))
    return(t(fold_change))
}

getAssay = function(se, colData_as_colname = 'Sample_Name', n=1){
    dat = unlist(assays(se,n))
    colnames(dat) = colData(se)[,colData_as_colname]
    return(dat)
}

reduce_genes_to_regulators = function(data, reg) {
    tfs = unique(reg$TF)
    tf_data <- data.frame(nrow=length(tfs),ncol=dim(data)[2])
    reg$fc_index = match(reg$name,rownames(data))
    # there are two tasks here that i am shoving into one and that is no good. first you want to get
    # the genes associated with a TF. and you may want to plot that. or stop there. THEN you want to
    # do some sort of agrigation of a TFs genes down into a single number to measure its activity.
}

rename_gene_names = function(matrix, t2g) {
  rownames(matrix) = t2g$name[match(rownames(matrix), t2g$target_id)]
  return(matrix)
}

average_duplicate_genes = function(data){
  dups = unique(rownames(data)[which(!isUnique(rownames(data)))])
  genes = unique(rownames(data))
  cond = colnames(data)
  ave_data = matrix(NA, nrow=length(genes), ncol = length(cond))
  ave_data = as.data.frame(ave_data)
  rownames(ave_data) = genes
  colnames(ave_data) = cond
  
  genes_not_in_dups = rownames(data)[which(isUnique(rownames(data)))]
  ave_data[genes_not_in_dups,] = data[genes_not_in_dups,]
  
  for(i in 1:length(genes)){
    if(genes[i] %in% dups){
      ave_data[genes[i],] = apply(data[which(rownames(data) == genes[i]),],2,function(x) mean(x, na.rm = T))
    }
  }
  return(ave_data)
}