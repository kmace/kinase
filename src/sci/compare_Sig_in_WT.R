# Data Import
source('../utils/load_libraries.R')
source('../utils/load_functions.R')
source('../utils/load_data.R')
# Get WT data
#raw_counts = raw_counts[ ,meta$Strain=='WT']
dds = dds[,colData(dds)$Strain == 'WT']
meta = meta[meta$Strain=='WT',]
dds = estimateSizeFactors(dds)
dds_counts = counts(dds, normalized=TRUE)
# Get Conditions available These might have drug or no drug
drugConditions = unique(meta[meta$Drug=="Cocktail",]$Condition)
noDrugConditions = unique(meta[meta$Drug=='None',]$Condition)

# Contrast genes
#We can perform this contrast in two ways. We can either see if the upregulated genes in Drug match those with no drug, or we can analyze each set individually, generating go terms for them.

#First thing is first. we need log foldchange with the appropriate baselines.

#Also here we are going to work on the mean for each condition. but I think it may be wise to look at each sample by sample comparrison individually. (use match for this)
getFC = function(counts, x, baseline){
 # base = rowMeans(counts[,meta$Condition==baseline],na.rm=T)
 # cond = rowMeans(counts[,meta$Condition==x],na.rm=T)
 base = apply(counts[,meta$Condition==baseline],1,function(x) median(x,na.rm=T))
 cond = apply(counts[,meta$Condition==x],1,function(x) median(x,na.rm=T))

 #bad = base==0 | cond==0
 #base = base[!bad]
 #cond = cond[!bad]
 #up = which(log2(cond) - log2(base) > 1)
 #down = which(log2(base) - log2(cond) > 1)
return(log2(cond) - log2(base))
}

contrast = c("Drug", "None", "Cocktail")
baseline = c("None_YPD_None", "None_YPD_Cocktail")

noDrug = do.call(cbind,lapply(noDrugConditions, function(x) getFC(dds_counts,x,baseline[1])))
drug = do.call(cbind,lapply(drugConditions, function(x) getFC(dds_counts,x,baseline[2])))

colnames(noDrug) = sapply(noDrugConditions,function(x) {a = strsplit(x, "_")[[1]][1]; if(a=="None"){a = strsplit(x, "_")[[1]][2]}; return(a)})
colnames(drug) = sapply(drugConditions,function(x) {a = strsplit(x, "_")[[1]][1]; if(a=="None"){a = strsplit(x, "_")[[1]][2]}; return(a)})

overlapping_conditions = colnames(noDrug)[colnames(noDrug) %in% colnames(drug)]
overlapping_conditions = overlapping_conditions[-which(overlapping_conditions=='YPD')]

plot_set = function(set) {
    set_up = set > 1
    set_up = set_up + 0
    upset(data.frame(set_up[complete.cases(set_up),]), nsets=dim(set_up)[2])
}

do_venn = function(cond) venn(list(drug=names(which(drug[,cond]>1)),no_drug = names(which(noDrug[,cond]>1))))
do_scatter = function(cond) plot(drug[,cond], noDrug[,cond])

dir.create('../../output/compare_Sig_in_WT.R')
pdf('../../output/compare_Sig_in_WT.R/venn.pdf')
lapply(overlapping_conditions, do_venn)
dev.off()
pdf('../../output/compare_Sig_in_WT.R/scatter.pdf')
lapply(overlapping_conditions, do_scatter)
dev.off()
# cond = "Fluconazole"; venn(list(drug=names(which(drug[,cond]>1)),no_drug = names(which(noDrug[,cond]>1))))
# cond = "Glucose Depletion"; plot(drug[,cond], noDrug[,cond])
