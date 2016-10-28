# ---- _3 ----
correlation_plot(deseq_counts)
correlation_plot(kallisto_counts)
correlation_plot(star_counts)


fc = filter_and_foldchange(deseq_counts)

correlation_plot(fc)

reg = load_regulation_matrix()
reg = add_reg_metadata(reg, t2g, fc)
plot_tf_expressions(expression, reg,file_name = 'test.pdf')
