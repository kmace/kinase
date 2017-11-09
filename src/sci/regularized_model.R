#load('intermediate/images/normalized_data.RData')
library(glmnet)
y = t(gene_expression_matrix)
X = as.data.frame(meta)[rownames(y),]

XX = model.matrix(~Strain + Condition + Strain:Condition, X)

fit = glmnet(x = XX, y=y, family='mgaussian')


