library(dplyr)
library(tidyr)
library(stringr)

files = read.table('files.txt', stringsAsFactors=F)
files = files$V1

all = do.call("rbind",
	lapply(files, function(fn) {
		if(length(readLines(fn)) > 12) {
			t = data.frame(Filename = fn,
										read.table(fn, skip=12));
			return(t)
			}
		}))

test = all %>%
	rename(TF = V7,
		left_pval = V13,
		right_pval = V14,
		two_tailed_pval = V15,
		corrected_left_pval = V20,
		corrected_right_pval = V21,
		corrected_two_tailed_pval = V22) %>%
	select(TF,
		left_pval,
		right_pval,
		two_tailed_pval,
		corrected_left_pval,
		corrected_right_pval,
		corrected_two_tailed_pval,
		Filename) %>%
	separate(Filename,
		c('dot', 'Kinase', 'Kinase_dir', 'f'),
		sep = '/') %>%
	select(-f, -dot) %>%
	separate(Kinase_dir, c('k', 's', 'Direction')) %>%
	select(-k, -s) %>%
	mutate(corrected_two_tailed_pval = as.numeric(str_sub(corrected_two_tailed_pval, 1, -2))) %>%
	select(Kinase, TF, corrected_two_tailed_pval) %>%
	arrange(corrected_two_tailed_pval)


test = test[test[,3]<0.01, ]
test = test %>% group_by(Kinase, TF) %>% summarize(pval = min(corrected_two_tailed_pval)) %>% ungroup

graph = tidyr::spread(test, TF, pval)
k = graph$Kinase
graph = as.data.frame(graph)
graph = graph[,-1]
rownames(graph) = k
bip = graph
bip = -log10(bip)
bip[is.na(bip)] = 0
bip = network(bip,
              matrix.type = "bipartite",
              ignore.eval = FALSE,
              names.eval = "weights")


col = c("actor" = "grey", "event" = "gold")
ggnet2(bip, color = "mode", palette = col, label = TRUE)
