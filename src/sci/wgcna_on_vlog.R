     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=4, 
     cex=.5)
colnames(matrix)
colnames(rlogm_all)
dist_all
dist_all[1]
dim(dist_all)
str(dist_all)
class(dist_all)
dist_all
dist_all@.S3Class
as.matrix(dist_all)
max(as.matrix(dist_all))
which(as.matrix(dist_all) == max(as.matrix(dist_all)))
dim(as.matrix(dist_all))
which(as.matrix(dist_all) == max(as.matrix(dist_all)),arr.ind=True)
which(as.matrix(dist_all) == max(as.matrix(dist_all)),arr.ind=TRUE)
apply(as.matrix(dist_all),1,mean)
which(as.matrix(dist_all) == max(as.matrix(dist_all)),arr.ind=TRUE)
max(apply(as.matrix(dist_all),1,mean))
t = apply(as.matrix(dist_all),1,mean))
t = apply(as.matrix(dist_all),1,mean)
which(t == max(t))
t = apply(as.matrix(dist_all),2,mean)
which(t == max(t))
dist_all[-29,-29]
as.matrix(dist_all)[-29,-29]
hclust(as.matrix(dist_all)[-29,-29])
hclust(as.dist(as.matrix(dist_all)[-29,-29]))
plot(hclust(as.dist(as.matrix(dist_all)[-29,-29])))
which(t %in% max(t))
which(t %in% sort(t)[1])
which(t %in% sort(t,decreasing=T)[1])
which(t %in% sort(t,decreasing=T)[1:3])
l = which(t %in% sort(t,decreasing=T)[1:3])
plot(hclust(as.dist(as.matrix(dist_all)[-l,-l])))
l = which(t %in% sort(t,decreasing=T)[1:4])
plot(hclust(as.dist(as.matrix(dist_all)[-l,-l])))
plot(hclust(as.dist(as.matrix(dist_all)[-l,-l])))
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=4, 
     cex=.5)
colnames(rlogm_all)
#colnames(rlogm_all) =
colnames(rld)
colData(rld)
colData(r)
ls()
colData(rld_all)
colData(rld_all)$Plate_Code
colnames(rlogm_all) = colData(rld_all)$Plate_Code
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=4, 
     cex=.5)
rlogm_all <- assay(rld)
lab = colData(rld_all)$Plate_Code
colnames(rlogm_all) = colData(rld_all)$Plate_Code
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
lab = lab[-leave]
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=4, 
     cex=.5)
lab
lab = colData(rld_all)$Media
colnames(rlogm_all) = colData(rld_all)$Plate_Code
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
lab = lab[-leave]
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=4, 
     cex=.5)
rlogm_all <- assay(rld)
lab = colData(rld_all)$Stress
colnames(rlogm_all) = colData(rld_all)$Plate_Code
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
lab = lab[-leave]
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=4, 
     cex=.5)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=2)
lab = colData(rld_all)$Drug
colnames(rlogm_all) = colData(rld_all)$Plate_Code
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
lab = lab[-leave]
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=2)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=2)
rlogm_all <- assay(rld)
lab = factor(colData(rld_all)$Plate_Code)
colnames(rlogm_all) = colData(rld_all)$Stress
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
lab = lab[-leave]
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=2)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=.5)
legend(2000,9.5, # places a legend at the appropriate place c(“Health”,”Defense”)
unique(meta$Plate_Code))
legend(2000,9.5, lab)
legend(lab)
legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8)
legend(1, 95, legend=unique(lab))
legend(10, 95, legend=unique(lab))
legend(10, 195, legend=unique(lab))
legend(10, 5, legend=unique(lab))
legend(0, 0, legend=unique(lab))
legend(10, 0, legend=unique(lab))
legend(100, 0, legend=unique(lab))
legend(100, 10, legend=unique(lab))
legend(100, 100, legend=unique(lab))
legend(100, 50, legend=unique(lab))
legend(100, 80, legend=unique(lab))
legend(100, 80, legend=unique(lab), color=hsv(as.numeric(lab)/max(as.numeric(lab))))
legend(100, 80, legend=unique(lab), col=hsv(as.numeric(lab)/max(as.numeric(lab))))
legend(100, 80, legend=unique(lab), fill=hsv(as.numeric(lab)/max(as.numeric(lab))))
legend('topright', legend=unique(lab), fill=hsv(as.numeric(lab)/max(as.numeric(lab))))
hsv(as.numeric(lab)/max(as.numeric(lab)))
hsv(as.numeric(unique(lab))/max(as.numeric(unique(lab))))
legend('topright', legend=unique(lab), hsv(as.numeric(unique(lab))/max(as.numeric(unique(lab)))))
legend('topright', legend=unique(lab), fill=hsv(as.numeric(unique(lab))/max(as.numeric(unique(lab)))))
pdf('distances.pdf')
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=.5)
legend('topright', 
       legend=unique(lab), 
       fill=hsv(as.numeric(unique(lab))/max(as.numeric(unique(lab)))))
dev.off()
pdf('distances.pdf', width=10, height=10)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=.5)
legend('topright', 
       legend=unique(lab), 
       fill=hsv(as.numeric(unique(lab))/max(as.numeric(unique(lab)))))
dev.off()
lab = factor(paste(colData(rld_all)$Stress,colData(rld_all)$Drug, colData(rld_all)$Media, sep='_' ))
lab
rlogm_all <- assay(rld)
rlogm_all <- assay(rld_all)
colnames(rlogm_all) = lab #colData(rld_all)$Stress
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
lab = lab[-leave]
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=.5)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=2)
plot(as.phylo(hc_subset), 
     type = "unrooted", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=.1)
     cex=.1)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
     label.offset=10, 
     cex=.5)
     type = "fan", 
history()
history(100)
rlogm_all <- assay(rld_all)
col_lab = factor(paste(colData(rld_all)$Stress,
                       colData(rld_all)$Drug, 
                       colData(rld_all)$Media, 
                       sep='_' ))
text_lab = factor(paste(colData(rld_all)$Stress,
                        colData(rld_all)$Drug, 
                        colData(rld_all)$Media,
                        colData(rld_all)$Strain,
                       sep='_' ))
#lab = factor(colData(rld_all)$Plate_Code)
colnames(rlogm_all) = text_lab #colData(rld_all)$Stress
dist_all <- dist(t(rlogm_all))
ave_dist = apply(as.matrix(dist_all),2,mean)
num_to_remove = 4
leave = which(t %in% sort(t,decreasing=T)[1:num_to_remove])
hc_subset = hclust(as.dist(as.matrix(dist_all)[-leave,-leave]))
col_lab = col_lab[-leave]
#dist_subset = dist(t(rlogm_all[,cutree(hc,h=150)==1]))
#hc_subset = hclust(dist_subset)
# plot(as.phylo(hc_subset), 
#      type = "unrooted", 
#      tip.color = hsv(as.numeric(lab)/max(as.numeric(lab))), 
#      label.offset=10, 
#      cex=.1)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(col_lab)/max(as.numeric(col_lab))), 
     label.offset=10, 
     cex=.5)
pdf('distances.pdf', width=10, height=15)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(col_lab)/max(as.numeric(col_lab))), 
     label.offset=10, 
     cex=.5)
dev.off()
pdf('distances.pdf', width=15, height=15)
plot(as.phylo(hc_subset), 
     type = "fan", 
     tip.color = hsv(as.numeric(col_lab)/max(as.numeric(col_lab))), 
     label.offset=10, 
     cex=.5)
dev.off()
dev.off()
q()
load('../input/images/normalized_data.RData')
library(dplyr)
meta$Condition
cbind(meta$Condition, meta$Strain_Code)
cbind(meta$Condition, factor(meta$Strain_Code))
meta %>% select(Condition, Strain_Code)
meta = meta %>% select(Condition, Strain_Code)
expr = vlog
plot(apply(vlog,1,mean), apply(vlog,1,sd))
apply(vlog,1,mean) > 6
expr = expr[apply(vlog,1,mean) > 6,]
model.matrix(meta)
model.matrix(~ Condition + Strain_Code, meta)
meta_matrix = model.matrix(~ Condition + Strain_Code, meta)
dim(meta_matrix)
dim(expr)
datExpr = t(expr)
datTraits = meta_matrix
library(WGCNA)
traitColors = numbers2colors(datTraits, signed = FALSE);
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
meta_matrix
colnames(meta_matrix)
enableWGCNAThreads()
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 10,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "kinaseTOM", 
                       verbose = 3)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                          signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = names(datTraits),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = greenWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1,1),
              main = paste("Module-trait relationships"))
names(datTraits)
datTraits
colnames(datTraits)
names(datTraits) = colnames(datTraits)
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = names(datTraits),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = greenWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1,1),
              main = paste("Module-trait relationships"))
?blockwiseModules
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "signed", minModuleSize = 10,
                       maxBlockSize = 1000, loadTOM = TRUE,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "kinaseTOM", 
                       verbose = 3)
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "signed", minModuleSize = 10,
                       maxBlockSize = 30000, loadTOM = TRUE,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "kinaseTOM", 
                       verbose = 3)
?blockwiseModules
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "signed", minModuleSize = 10,
                       maxBlockSize = 30000, loadTOM = TRUE,
                       deepSplit = 4, 
                       reassignThreshold = 0,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "kinaseTOM", 
                       verbose = 3)
table(net$colors)
meta_matrix = model.matrix(~ Condition + Strain_Code + 0, meta)
meta_matrix
colnames(meta_matrix)
?model.matrix
meta$Condition
colnames(meta_matrix)
lapply(meta,length)
lapply(meta,function(x)length(levels(x)))
lapply(meta,function(x)length(unique(x)))
dim(meta_matrix)
meta_matrix = model.matrix(~ Condition + Strain_Code + 0, meta, contrast.arg=lapply(meta,contrasts, contrasts=FALSE))
dim(meta_matrix)
library(caret)
install.packages('caret')
library(caret)
dummyVars(meta)
dummyVars(~., meta)
meta
dummyVars(~., data=meta, levelsOnly=T)
test = dummyVars(~., data=meta, levelsOnly=T)
test
as.matrix(test)
predict(test, meta)
test = predict(test, meta)
dim(test)
colnames(test)
class(meta$Condition)
class(meta$Strain_Code)
meta$Strain_Code = factor(meta$Strain_Code)
class(meta$Strain_Code)
meta
test = predict(dummyVars(~., meta), meta)
test
dim(test)
colnames(test)
test = predict(dummyVars(~., meta, levelsOnly=T), meta)
test
colnames(test)
meta_matrix = predict(dummyVars(~., meta, levelsOnly=T), meta)
names(meta_matrix)
names(meta_matrix) = colnames(meta_matrix)
datTraits
dim(datTraits)
datTraits
datTraits = meta_matrix
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                          signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = names(datTraits),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = greenWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1,1),
              main = paste("Module-trait relationships"))
pdf('test.pdf')
labeledHeatmap(Matrix = moduleTraitCor,
              xLabels = names(datTraits),
              yLabels = names(MEs),
              ySymbols = names(MEs),
              colorLabels = FALSE,
              colors = greenWhiteRed(50),
              textMatrix = textMatrix,
              setStdMargins = FALSE,
              cex.text = 0.5,
              zlim = c(-1,1),
              main = paste("Module-trait relationships"))
dev.off()
ls()
?history
savehistory(file='wgcna_on_vlog.Rhistory')
