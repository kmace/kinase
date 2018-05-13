writeGCT <- function (data,
                      row_metadata,
                      col_metadata,
                      filename) {
    nrows = dim(data)[1]
    ncols = dim(data)[2]
    row_metadata_fields = dim(row_metadata)[2]
    col_metadata_fields = dim(col_metadata)[2]
    f <- file (filename, "w")
    on.exit (close (f))
    cat ("#1.3", "\n", file = f, append = TRUE, sep = "")
    cat (nrows, "\t",
         ncols, "\t",
         row_metadata_fields, "\t",
         col_metadata_fields, "\n", file = f,
         append = TRUE, sep = "")
    # First row
    cat ("Gene/id", file = f, append = TRUE, sep = "")
    names <- colnames(row_metadata)
    for (j in 1:length(names)) {
        cat ("\t", names[j], file = f, append = TRUE, sep = "")
    }
    names <- colnames(data)
    for (j in 1:length (names)) {
        cat ("\t", names[j], file = f, append = TRUE, sep = "")
    }
    cat ("\n", file = f, append = TRUE, sep = "")
    # Build col metadata
    names = colnames(col_metadata)
    nas = matrix(NA, nrow= col_metadata_fields, ncol=row_metadata_fields  )
    col_metadata_matrix = cbind(names, nas, t(as.data.frame(lapply(col_metadata,as.character))))
    write.table(col_metadata_matrix, file = f, append = TRUE, quote = FALSE, sep = "\t",
                 eol = "\n", col.names = FALSE, row.names = FALSE)
    # Build row metadata
    names = rownames(data)
    m <- cbind (names, row_metadata, data)
    write.table (m, file = f, append = TRUE, quote = FALSE, sep = "\t",
                 eol = "\n", col.names = FALSE, row.names = FALSE)
}
