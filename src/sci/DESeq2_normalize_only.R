#! /usr/bin/Rscript

###
###  Run DESeq2 from the command line
###  George Bell - Bioinformatics and Research Computing, Whitehead Institute
###
###  USAGE: DESeq2_normalize_only.R inputCounts OutputFile 
###  EX: DESeq2_normalize_only.R InputDEseq2.txt DESeq2_output.txt
###
###
###  12 August 2016 [GWB]:
###    Modify DESeq2 code for Rscript syntax 
###      and make PCA figures with similar code as DESeq_normalize_only.R
###

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./DESeq2_normalize_only.R" }

if (length(commandArgs()) < (7 - offset))
{
	message("\nNormalize a matrix of raw counts with DEseq2.")
	message(paste("USAGE:   ", this.script, "inputCounts OutputFile"))
	message(paste("Example: ", this.script, "InputDEseq2.txt DESeq2_output.2016.txt\n"))
	q()
}

# Load library without messages
message("\nLoading DESeq2 and other required R packages ....")
suppressMessages(library(DESeq2))
suppressMessages(library(genefilter))
suppressMessages(library(RColorBrewer))
suppressMessages(library(lattice))
suppressMessages(library(ggrepel))

# First argument is input file
input.filename = commandArgs()[6 - offset]
# Second argument is output file
output.filename = commandArgs()[7 - offset]

message(paste("Input filename is", input.filename))

# Load columns of expression counts
counts = read.delim(input.filename, row.names=1)

# pca.figure.filename = paste(input.filename, "PCA.pdf", sep=".")
# paste("PDF file of PCA figure is", pca.figure.filename)


# Make a DESeqDataSet
samples = factor(colnames(counts))
# Use a design stating that there are no replicates
# If there is replication, replace 'design' with the real design
message("Treating all samples as independent (so assuming no replication)....\n")
dds = DESeqDataSetFromMatrix(countData=counts, colData=DataFrame(samples), design=~1)

###  MAJOR COMMAND:
# Estimate size factors and dispersions and fit generalized linear model
# For designs without replication, fitType has no influence on size factors
# dds = DESeq(dds, fitType="local")
dds = DESeq(dds, fitType="mean")

# Print size factors
message("\nSize factors:")
sizeFactors(dds)
message()

# Add normalized counts for this version (GB - 19 Mar 2013)
counts.normalized = round(t(t(counts(dds))/sizeFactors(dds)), 2)
colnames(counts.normalized) = paste(colnames(counts), "norm", sep=".")

# Print output (including norm counts)
output.table = cbind(rownames(dds), counts, counts.normalized)
colnames(output.table)[1] = "Feature.ID"
write.table(output.table, file=output.filename, sep="\t", quote=F, row.names=F)


# We may need to modify the function plotPCA() because it makes a max of 12 colors (so >12 samples cause problems)
plotPCA_BaRC_DESeq2 = function (x, intgroup = "condition", ntop = 500, aspect="iso")
{
	# Need to convert SummarizedExperiment into a plain old matrix
	x = assay(x)
    rv = rowVars(x)
    groups = colnames(x)
    select = order(rv, decreasing = TRUE)[seq_len(ntop)]
    pca = prcomp(t(x[select, ]))
	colours = brewer.pal(ncol(x), "Paired")
    xyplot(PC2 ~ PC1, groups = groups, data = as.data.frame(pca$x), 
        pch = 16, cex = 2, aspect = aspect, col = colours, main = draw.key(key = list(rect = list(col = colours), 
        text = list(groups), rep = FALSE)))
}

ggplot.pca = function (x, y, fac, pca.x.name, pca.y.name, ntop=500)
{
	# Set up colors and shapes of points
    if (length(fac) >= 3)
        colours = rainbow(length(fac))
	ggplot() +
	geom_point(aes(x, y), color = colours, size=5) + 
	labs(x = pca.x.name) +
	labs(y = pca.y.name) + 
	theme_classic(base_size = 16) +
	ggtitle(paste("Using top", ntop, "features")) + 
	geom_text_repel(aes(x, y, label = fac)) + 
	theme_bw()	# To add axis lines
	# coord_equal(ratio = 1)	# Make isometric axes (so each scale is represented by the same distance)
	# geom_text(aes(x, y), label=fac, hjust=0, nudge=nudge, check_overlap = TRUE)
}

plotPCA_BaRC_v3 = function (x, ntop = 500)
{
	# Need to convert SummarizedExperiment into a plain old matrix
	x = assay(x)
    rv = rowVars(x)
	rv = genefilter::rowVars(x)
	select = order(rv, decreasing = TRUE)[seq_len(ntop)]
	pca = prcomp(t(x[select, ]))
	
	#add this line to get the weights of the genes on each component
	write.table(pca$rotation, file=paste(input.filename, "pca_genes.using.top", ntop, "txt", sep="."),sep="\t", quote=F)
	# Get the weights of the samples on each component
	write.table(pca$x, file=paste(input.filename, "pca_samples.using.top", ntop, "txt", sep="."),sep="\t", quote=F)

	# Eigenvalues
	eig = (pca$sdev)^2
	variances = round(eig*100/sum(eig), 2)
	names(variances) = colnames(pca$x)

	par(las=1)
	print(ggplot.pca(pca$x[,"PC1"], pca$x[,"PC2"], rownames(pca$x), paste("PC1:   ", variances[1], "% variance", sep=""), paste("PC2:   ", variances[2], "% variance", sep=""), ntop))
	print(ggplot.pca(pca$x[,"PC1"], pca$x[,"PC3"], rownames(pca$x), paste("PC1:   ", variances[1], "% variance", sep=""), paste("PC3:   ", variances[3], "% variance", sep=""), ntop))
	print(ggplot.pca(pca$x[,"PC2"], pca$x[,"PC3"], rownames(pca$x), paste("PC2:   ", variances[2], "% variance", sep=""), paste("PC3:   ", variances[3], "% variance", sep=""), ntop))

	message(paste("Variances (as percentage) for top", ntop, "features"))
	print(variances)
	# barplot(eig, ylab="eigenvalues", main=paste("Eigenvalues for top", ntop, "features"))
	barplot(variances, ylab="variance (as percentage)", main=paste("Variances (as percentage) for top", ntop, "features"))
}


### Do VST transformation (required to  for PCA plot to remove the dependence of the variance on the mean)
# From manual:
# blind=TRUE should be used for comparing samples in an manner unbiased by prior information on samples, for example to perform sample QA (quality assurance). 
# blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made.
dds.vst = varianceStabilizingTransformation(dds, blind=TRUE, fitType="mean")
# rld = rlogTransformation(dds, blind=TRUE, fitType="mean")

### Run Principal Components Analysis
pca.figure.filename = paste(input.filename, "PCA_etc.pdf", sep=".")
pdf(pca.figure.filename, w=8.5, h=8.5)
# Fill the page
# plotPCA_BaRC_DESeq2(dds.vst, aspect="fill")
# Make both axes on the same scale
# plotPCA_BaRC_DESeq2(dds.vst, aspect="iso")

ntops = c(300, 500, 1000)
for (i in 1:length(ntops))
{
	plotPCA_BaRC_v3(dds.vst, ntop=ntops[i])
}

num.samples = length(samples)
par(las=2, mai=c(2,1,1,1))
boxplot(log(cbind(counts, counts.normalized) + 1, 2), ylab="log2 ( counts + 1 )", col=c(rep("gray",num.samples),rep("wheat",num.samples)))

# Close output file
foo = dev.off

message(paste("\nPCA plots are in", pca.figure.filename))
message(paste("All done -- see", output.filename, "for normalized counts.\n"))


###################

