load('../../intermediate/images/normalized_data.RData')
library(edgeR)
library(stringr)
library(dplyr)

dge = readDGE(str_replace(meta$star_path, 'input','intermediate' ), columns = c(1,3), labels = meta$Experiment)
dge = dge[-c(1:3),]
rownames(dge$counts) = t2g$name[match(rownames(dge$counts),t2g$Gene)]
genes = t2g[match(rownames(dge),t2g$name),]
x = dge
x$genes <- genes
#x
#cpm <- cpm(x)
#lcpm <- cpm(x, log=TRUE)

#keep.exprs <- rowSums(lcpm>1)>=50
#x <- x[keep.exprs,, keep.lib.sizes=FALSE]
#dim(x)
library(RColorBrewer)
# nsamples <- ncol(x)
# col <- brewer.pal(nsamples, "Paired")
# par(mfrow=c(1,2))
# plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
# main="", xlab="")
# title(main="A. Raw data", xlab="Log-cpm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
# den <- density(lcpm[,i])
# lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", samplenames, text.col=col, bty="n")
# lcpm <- cpm(x, log=TRUE)
# plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
# main="", xlab="")
# title(main="B. Filtered data", xlab="Log-cpm")
# abline(v=0, lty=3)
# for (i in 2:nsamples){
# den <- density(lcpm[,i])
# lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", samplenames, text.col=col, bty="n")
x <- calcNormFactors(x, method = "TMM")
x$samples = left_join(x$samples %>% mutate(Experiment = rownames(x$samples)), meta)
rownames(x$samples) = x$samples$Experiment
glMDSPlot(lcpm, labels=x$samples$Strain,
groups=x$samples %>% select(Strain, Condition), launch=T)

meta$Strain %<>% droplevels()
meta$Condition %<>% droplevels()


design = model.matrix(~ 0 + Condition + 0 + Strain + Strain:Condition, meta)

#all.zero <- apply(design, 2, function(x) all(x==0))
#idx <- which(all.zero)
#design <- design[,-idx]


colnames(design) = colnames(design) %>% str_replace('Strain','') %>% str_replace('Condition','')  %>% str_replace(':','_') %>%  str_replace(' ','')

makeContrasts(Salt,
              Salt_PBS2 - Salt,
              Salt_HOG1 - Salt,
              levels=colnames(design)) -> contr.matrix

makeContrasts(Tunicamycin,
IRE1,
Tunicamycin_IRE1,
Tunicamycin_IRE1 - Tunicamycin,
Tunicamycin_IRE1 - (Tunicamycin + IRE1),
levels=colnames(design)) -> contr.matrix


v <- voom(x, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

efit <- eBayes(vfit)
dt = decideTests(efit)

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


de.common <- which(dt[,1]!=0 & dt[,2]!=0 & dt[,3]!=0)
de.common
vennDiagram(dt, circle.col=c("turquoise", "salmon", "gold"))
vfit
tfit <- treat(vfit, lfc=.5)
dt <- decideTests(tfit)
summary(dt)
write.fit(tfit, dt, file="results.txt")
getwd()
topTreat(tfit, coef=1, n=Inf)
topTreat(tfit, coef=3, n=Inf)
basal.vs.lp <- topTreat(tfit, coef=3, n=Inf)
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],
xlim=c(-8,13))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],
xlim=c(-8,13))
plotMD(tfit, column=3, status=dt[,3], main=colnames(tfit)[3],
xlim=c(-8,13))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1],
xlim=c(-8,13))
basal.vs.lp
basal.vs.lp %>% arrange(desci(logFC))
basal.vs.lp %>% arrange(desc(logFC))
basal.vs.lp %>% rownames_to_column('Gene') %>% arrange(desc(logFC)) %>%
left_join(t2g)
makeContrasts(Salt - YPD, PBS2, Salt_PBS2 - Salt, Salt_PBS2 - (Salt + PBS2), Salt_PBS2 - PBS2, levels=colnames(design)) -> contr.matrix
makeContrasts(Salt - YPD, PBS2, Salt_PBS2 - (Salt + PBS2), (Salt_PBS2 + Salt + PBS2) - (Salt + PBS2), Salt_PBS2 - PBS2, levels=colnames(design)) -> contr.matrix
contr.matrix
rowSums(abs(contr.matrix))
rowSums(abs(contr.matrix))!=0
contr.matrix[rowSums(abs(contr.matrix))!=0,]
#vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
nrow(contr.matrix)
vfit$coefficients %>% dim()
vfit$coefficients = NULL
#vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean−variance trend")
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)
head(dt)
dt[,c(1,3,4)]
vennDiagram(dt[,c(1,3,4)])
dt %>% rownames_to_column('Gene') %>% left_join(t2g)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(.[[4]] != 0)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(.[[3]] != 0)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(.[[4]] != 0)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(.[[5]] != 0)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(.[[4]] != 0)
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(name == 'FUS1')
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(name == 'AGA1')
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(name == 'AGA2')
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(name == 'PUA4')
dt %>% as.data.frame() %>% rownames_to_column('Gene') %>% left_join(t2g) %>% filter(name == 'PAU4')
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
counts=x$counts, groups=x$samples %>% select(Condition, Strain), launch=T)
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
counts=x$counts, groups=x$samples %>% pull(Condition), launch=T)
glMDPlot(tfit, coef=4, status=dt, main=colnames(tfit)[4],
counts=x$counts, groups=x$samples %>% pull(Condition), launch=T)
glMDPlot(tfit, coef=4, status=dt, main=colnames(tfit)[4],
counts=x$counts, groups=x$samples %>% pull(Condition), launch=T, side.main = "name")
glMDPlot(tfit, coef=4, status=dt, main=colnames(tfit)[4],
counts=x$counts, groups=x$samples %>% select(Condition, Strain), launch=T, side.main = "name")
?glMDPlot
?camera
camera(v, which(dt[,4] == 1, design = design)
)
camera(v, which(dt[,4] == 1), design = design)
cameraPR(efit$t[,2], list(set1=index1,set2=index2))
efit$t
cameraPR(efit$t[,2], list(set1=c("YER026C", "YER027C")))
modules
cameraPR(efit$t[,2], list(set1=c("YER026C", "YER027C")))
cameraPR(efit$t[,3], list(set1=c("YER026C", "YER027C")))
cameraPR(efit$t[,4], list(set1=c("YER026C", "YER027C")))
