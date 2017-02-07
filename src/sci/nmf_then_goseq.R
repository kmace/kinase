library(NMF)
res = nmf(val,12)
fit(res)
val.hat = fitted(res)
dim(val.hat)
dim(val)
hist(val - val.hat)
summary(res)
# get matrix W
w <- basis(res)
dim(w)
## [1] 200 3
# get matrix H
h <- coef(res)
dim(h)
image(w)
image(h)
s <- extractFeatures(res)

get_go_terms_from_list = function(query_list, background_list) {
    library(goseq)
    genes = numeric(length(background_list))
    names(genes) = background_list
    genes[query_list] = 1
    pwf = nullp(genes,"sacCer2","ensGene")
    GO.wall=goseq(pwf,gene2cat=as.list(org.Sc.sgdGO2ALLORFS))
    return(GO.wall)
}

get_go_terms_from_list(s[[11]]) %>% head
