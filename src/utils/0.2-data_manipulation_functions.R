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

reduce_genes_to_regulators = function(data, reg) {
    tfs = unique(reg$TF)
    tf_data <- data.frame(nrow=length(tfs),ncol=dim(data)[2])
    reg$fc_index = match(reg$name,rownames(data))
    # there are two tasks here that i am shoving into one and that is no good. first you want to get
    # the genes associated with a TF. and you may want to plot that. or stop there. THEN you want to
    # do some sort of agrigation of a TFs genes down into a single number to measure its activity.
}
