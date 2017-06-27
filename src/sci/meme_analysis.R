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
	       num_genes = V10,
		left_pval = V13,
		right_pval = V14,
		two_tailed_pval = V15,
		corrected_left_pval = V20,
		corrected_right_pval = V21,
		corrected_two_tailed_pval = V22) %>%
	select(TF,
	       num_genes,
		left_pval,
		right_pval,
		two_tailed_pval,
		corrected_left_pval,
		corrected_right_pval,
		corrected_two_tailed_pval,
		Filename) %>%
	separate(Filename,
		c('dot', 'exp', 'f'),
		sep = '/') %>%
	select(-f, -dot) %>%
	separate(exp, c('Kinase', 'pval', 'method', 'param', 'junk', 'Direction'), sep='_') %>%
	select(-junk) %>%
	mutate(corrected_two_tailed_pval = as.numeric(str_sub(corrected_two_tailed_pval, 1, -2))) %>%
	#select(Kinase, TF, corrected_two_tailed_pval) %>%
	arrange(corrected_two_tailed_pval)


test %>% filter(TF=='MSN2') %>% ggplot(aes(x = Kinase, y = -log10(corrected_two_tailed_pval))) + geom_point(aes(color = Direction))
test %>% filter(TF=='RLM1') %>% ggplot(aes(x = Kinase, y = -log10(corrected_two_tailed_pval))) + geom_point(aes(color = Direction))
test %>% filter(TF=='HAC1') %>% ggplot(aes(x = Kinase, y = -log10(corrected_two_tailed_pval))) + geom_point(aes(color = Direction))
test %>% filter(TF=='PHO4') %>% ggplot(aes(x = Kinase, y = -log10(corrected_two_tailed_pval))) + geom_point(aes(color = Direction))

test %>% ggplot(aes(x = Kinase, y = -log10(corrected_two_tailed_pval))) + geom_point() + facet_wrap(~TF)

test %>% group_by(TF, Kinase) %>% summarise(count = n()) %>% ggplot(aes(x=TF, y=Kinase)) + geom_tile(aes(fill = count))

do_plot = function(tf, kinase){
  test %>% 
    filter(TF==tf, Kinase == kinase) %>% 
    ggplot(aes(x = pval, y = param, fill=-log10(corrected_two_tailed_pval))) + 
    geom_tile() + facet_grid(~method) + ggtitle(paste(kinase, 'interacts with', tf))
}

do_plot('MSN2', 'TPK123')

test %>% group_by(TF, Kinase) %>% 
  summarise(count = n()) %>% 
  #filter(count > 15) %>% 
  ggplot(aes(x=TF, y=Kinase)) + 
    geom_tile(aes(fill = count)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=5))

net_data = test %>% filter(param == 200, pval == 0.01, corrected_two_tailed_pval < 0.01) %>% select(TF, Kinase, corrected_two_tailed_pval)

test = test[test[,3]<0.01, ]
test = test %>% group_by(Kinase, TF) %>% summarize(pval = min(corrected_two_tailed_pval)) %>% ungroup

graph = tidyr::spread(net_data, TF, corrected_two_tailed_pval)
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
