load('../../input/images/normalized_data.RData')
load('../../input/images/dictionary.RData')

library(DESeq2)
target = assay(vsd)

dictionary = dictionary[rownames(dictionary) %in% rownames(target),]
target = target[rownames(target) %in% rownames(dictionary),]
gm = apply(target,1,mean)
target = apply(target,2,function(x) x - gm)

dictionary = dictionary[match(rownames(target) , rownames(dictionary)),]

library(glmnet)
fit1 = glmnet(dictionary,target[,1])

get_10_c = function(y) {
  c = coef(glmnet(dictionary,y,dfmax=10))
  last_c = dim(c)[2]
  return(c[,last_c])
}

all_c = apply(target,2,get_10_c)
top_features = names(sort(apply(all_c,1,function(x) sum(abs(x))),decreasing=T)[1:50])
b = all_c[top_features,]
d = dictionary[,rownames(b)]
m = mldivide(d,target)
d3heatmap(m,col=viridis(20))
