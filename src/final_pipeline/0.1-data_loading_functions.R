# ---- _0.1 ----
load_sample_meta_data = function() {
	meta = read.csv('../../input/meta/Sample_Meta.csv', header=T, stringsAsFactors=FALSE)
	meta =  meta %>%
					mutate(Strain = relevel(factor(Strain),'WT'),
								 Drug = relevel(factor(Drug),'None'),
								 Stress = relevel(factor(Stress),'None'),
								 Media = relevel(factor(Media),'YPD'),
								 Experimenter = relevel(factor(Experimenter),'Kieran')) %>%
					unite(Experiment, c(Stress, Media, Drug, Strain_Code), remove = FALSE) %>%
					mutate(Experiment = relevel(factor(Experiment), 'None_YPD_None_WT1'))
		return(meta)
}

load_seq_meta_data = function() {
	# lane_annotation files look differnet, in terms of the number of feilds. this will take some thinking to figure out, need to for loop
	# over the annotation files, and pull in relivant fields, or maybe we just exclude previous data.
	meta_files = list.files('../../input/meta/all_seq_meta/',pattern = '*lane*')
	meta <- do.call("rbind.fill",
					lapply(meta_files,
						   function(fn) {
								 t = data.frame(Filename=fn,
							   					   		read.table(file.path('../../input/meta/all_seq_meta',fn),
												   			header=T,
												   			sep='\t',
												   			stringsAsFactors=FALSE))
								#t = t[,-which(colnames(t) %in% c("NumberOfPeaks","X.RiP","Strain", "Drug","Condition","concbyQubit.nM..DF."))]
						   	return(t)
							}
						  )
					)
	meta = meta %>%
		mutate(tophat_path = file.path('../../intermediate/counts/tophat', paste(SampleName, '_gene_counts.txt',sep='')),
					 star_path = file.path('../../intermediate/counts/star', paste(SampleName, '_sequence.txt.gz_ReadsPerGene.out.tab',sep='')),
					 star_column = if_else(Adaptor=="TruSeq-RNA",3,2),
					 kallisto_path = file.path('../../intermediate/kallisto/kallisto_output', paste(SampleName,'_sequence.txt.gz',sep='')) ) %>%
		rename(sample = SampleName)
		return(meta)
}

get_sample_meta = function() {
		sample_meta = load_sample_meta_data()
		seq_meta = load_seq_meta_data()
		meta = right_join(sample_meta, seq_meta, by = c('Sample_Name' = 'sample'))
		return(meta)
}

load_transcripts_to_genes = function() {
		mart <- biomaRt::useMart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
		t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
											 "ensembl_gene_id",
											 "external_gene_name",
											 "description"), mart = mart)

		t2g <- dplyr::rename(t2g,
			target_id = ensembl_transcript_id,
			ens_gene = ensembl_gene_id,
			ext_gene = external_gene_name)
		t2g <- t2g %>% dplyr::mutate(name=ifelse(ext_gene=='',ens_gene,ext_gene))
		return(t2g)
}

load_count_file = function(path, type = 'tophat', star_col=NA) {
		counts_table = read.table(path, row.names=1)
		d = dim(counts_table)
		if(type == 'star') {
		    if(is.na(star_col)){
		      stop("if you want to load reads from star, you need to indicate what column to use")
		    }
				num_non_genes = 4
				counts_table = counts_table[-c(1:num_non_genes), star_col, drop=F]
				# non_gene_location = 'begining'
		} else if(type == 'tophat') {
				num_non_genes = 5
				counts_table = counts_table[-c((d[1]-num_non_genes+1):d[1]), , drop=F]
				# non_gene_location = 'end'
		} else {
				stop('Bad file type argument')
		}
		return(counts_table)
}

load_count_matrix = function(meta, type = 'tophat'){
		# The meta file must contain at least the seq meta data.
		paths = switch(type, tophat = meta$tophat_path, star = meta$star_path)
		col = switch(type, tophat = NA, star = meta$star_column)
		count_matrix <- do.call("cbind",
										mapply(function(p, c) load_count_file(p, type, c),
										       paths,
										       col, SIMPLIFY = FALSE))
		colnames(count_matrix) = meta$Sample_Name
		return(count_matrix)
}

load_gene_list_from_file = function(file) {
		return(read.table(file, header=F, stringsAsFactors=F)[,1])
}

load_binary_yeastract_regulation_table = function(file = '../data/yeastract/RegulationTwoColumnTable_Documented_2013927.tsv', sep = ';') {
		reg = fread(file,
								header=T,
								sep=sep)
		# regulation = fread('../data/yeastract/RegulationMatrix_Documented_2013927.csv',
		#												 row.names=1,
		#												 quote='',
		#												 header=T,
		#												 sep=';')
		colnames(reg) = c('TF', 'Target')
		return (reg)
}

load_probBinding_yetfasco_regulation_table = function(file = '../../input/reference/yetfasco/20120129_allMotifData1.02.rdat', use_expert=TURE) {
		load(file)
		library(rvest)
		library(stringr)
		library(tidyr)
		txt = paste('<table>', str_c(datasetNames, collapse=' '), '</table>')
		html = read_html(txt)
		exp_table = html_table(html)[[1]]
		regT = data.frame(t(dataMat))
		regT[,'expert'] = FALSE
		regT[expert,'expert'] = TRUE
		regT[,'TF'] = exp_table[,2]
		regT[,'exp_name'] = exp_table[,1]
		reg = gather(regT, Target, log_prob_bind, -expert, -TF, -exp_name)
		return(reg)
}

load_pval_rickYong_regulation_table = function(file = '~/Desktop/Datasets/106_pvalbygene_ypd_v9.0.csv') {
		regT = read.csv(file,header=T, stringsAsFactors=F)
		reg = gather(regT,TF,pval, -ens_name,-ext_name,-description) %>%
		rename(Target=ext_name)
		reg$Target[reg$Target=="#REF!"] = reg$ens_name[reg$Target=="#REF!"]
		return(reg)
}

load_slueth_object = function(model, meta, t2g) {
		so <- sleuth_prep(meta, model, target_mapping = t2g, aggregation_column = 'ens_gene')
		so <- sleuth_fit(so)
		return(so)
}

load_deseq_object = function(model, meta, gene_count_matrix) {
		dds<-	DESeqDataSetFromMatrix(countData= gene_count_matrix, colData= meta, model) #~ Strain + Drug + Condition)
		dds <- DESeq(dds)
		dds <- estimateSizeFactors(dds)
		return(dds)
}

load_gasch = function() {
  gasch = read.table('../../input/reference/gasch_complete_dataset.txt', header=T, sep='\t', quote="", row.names=1, stringsAsFactors=F)
  #g_info = gasch[,"NAME"]
  gmeta = read.csv('../../input/meta/gasch_meta.csv')
  Stress_after = as.character(gmeta$Stress)
  Stress_after[gmeta$Time<10] = 'Pre-response'
  gmeta$Stress_after = factor(Stress_after)

  gx = gasch[,-c(1,2)]
  gx = impute.knn(as.matrix(gx))$data
  gx = gx - rowMeans(gx)
  gasch = list()
  gasch$data = gx
  gasch$meta = gmeta
  return(gasch)
}

load_megaYeast = function() {
  mega = read.table('../../input/reference/Mega_YeastX.txt',
                     header=T,
                     sep='\t',
                     quote="",
                     comment.char = '',
                     row.names=1,
                     stringsAsFactors=F)
  #g_info = gasch[,"NAME"]
  #mmeta = read.csv('../../input/meta/gasch_meta.csv')
  #Stress_after = as.character(gmeta$Stress)
  #Stress_after[gmeta$Time<10] = 'Pre-response'
  #gmeta$Stress_after = factor(Stress_after)

  mx = mega[,-1]
  num_missing_genes = apply(mx,2,function(x) sum(is.na(x)))
  mx = mx[,num_missing_genes<3000]
  num_missing_genes = apply(mx,1,function(x) sum(is.na(x)))
  mx = mx[num_missing_genes<400, ]

  mx = impute.knn(as.matrix(mx))$data
  mx = mx - rowMeans(mx)
  mega = list()
  mega$data = mx
  #mega$meta = mmeta
  return(mega)
}

load_pronk = function() {
  gset <- getGEO("GSE11452", GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset) > 1) idx <- grep("GPL90", attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]
  pmeta = pData(gset)
  pmeta$Condition = as.factor(substr(as.character(pmeta$title), start=1, stop=30))
  px = exprs(gset)
  p2t = read.table('../../input/reference/GPL_to_ORF.tsv', fill=T, header=T, stringsAsFactors=F)
  rownames(px) = p2t[,2]
  px = px[rownames(px) != '',]
  px = log2(px + 0.0001) # smallest nonzero is 2^-12
  px = impute.knn(px)$data
  px = px - rowMeans(px)
  px = average_duplicate_genes(px)

  # rv = genefilter::rowVars(gx)
  # select = order(rv, decreasing = TRUE)[seq_len(4000)]
  # gx = gx[select, ]
  colnames(px) = pmeta$title
  pronk = list()
  pronk$data = px
  pronk$meta = pmeta
  return(pronk)
}
