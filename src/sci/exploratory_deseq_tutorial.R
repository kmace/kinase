dds = deseq_object
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)
rld <- rlog(dds, blind=FALSE)

#Scatterplot of transformed counts from two samples. Shown are scatterplots using the log2 transform of normalized counts (left side) and using the rlog (right side).

par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)



sampleDists <- dist( t( assay(rld) ) )



#Heatmap of sample-to-sample distances using the rlog-transformed values.
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Poisson Distance (Witten 2011), implemented in the PoiClaClu package. This measure of dissimilarity between counts also takes the inherent variance structure of counts into consideration when calculating the distances between samples
poisd <- PoissonDistance(t(counts(dds)))

 samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- rld$sample
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)


plotPCA(rld, intgroup = c("Condition", "kinase_deactivated"))
plotPCA(rld, intgroup = c("Condition", "Strain","Drug"))
plotPCA(rld, intgroup = c("Condition"))


mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=kinase_deactivated,shape=Condition)) + geom_point(size=3) +
coord_fixed()


dds <- DESeq(dds)
res = results(dds)
summary(res)
table(res$padj < 0.1)
