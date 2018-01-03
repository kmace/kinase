library(RColorBrewer)
library(grDevices)
library(ggplot2)


kinase_colors = c(
'#ffffff',
'#003B4A',
'#1f77b4',
'#aec7e8',
'#ff7f0e',
'#ffbb78',
'#637939',
'#38C697',
'#2ca02c',
'#98df8a',
'#FB4124',
'#830303',
'#d62728',
'#ff9896',
'#393b79',
'#9467bd',
'#c5b0d5',
'#3f1818',
'#8c564b',
'#c49c94',
'#e377c2',
'#f7b6d2',
'#7f7f7f',
'#c7c7c7',
'#000000',
'#bcbd22',
'#dbdb8d',
'#17becf',
'#9edae5')

names(kinase_colors) = c(
"WT",
"ATG1",
"CDC15",
"CDC5",
"CDC7",
"CTK1",
"FUS3",
"HOG1",
"IRE1",
"KIN1",
"KSP1",
"KSS1",
"MRK1",
"PBS2",
"PHO85",
"PKC1",
"RIM11",
"RIM15",
"SCH9",
"SLT2",
"SNF1",
"SSN3",
"STE11",
"TOR2",
"TPK123",
"YAK1",
"YGK3",
"YPK1",
"YPK3")

condition_colors = RColorBrewer::brewer.pal(10,"Set3")
names(condition_colors) = sort(unique(module_expression$Condition))

condition_color_scale = scale_color_manual(name = 'Condition', values = condition_colors)
kinase_color_scale = scale_color_manual(name = 'Kinase', values = kinase_colors)
strain_color_scale = scale_color_manual(name = 'Strain', values = kinase_colors)

condition_fill_scale = scale_fill_manual(name = 'Condition', values = condition_colors)
kinase_fill_scale = scale_fill_manual(name = 'Kinase', values = kinase_colors)
strain_fill_scale = scale_fill_manual(name = 'Strain', values = kinase_colors)

master_col = list(Condition = condition_colors,
                  Kinase = kinase_colors)

coolwarm_hcl <- colorspace::diverge_hcl(11,
                                        h = c(250, 10),
                                        c = 100,
                                        l = c(37, 88),
                                        power = c(0.7, 1.7))[-c(1, 11)]



theme_Publication <- function(base_size=14, base_family='') {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))

}
