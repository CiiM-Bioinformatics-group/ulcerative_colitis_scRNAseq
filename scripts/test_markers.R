#setwd("/path/to/projects")
options(future.globals.maxSize = 10737418240)  # 10240*1024^2 = 10737418240  for 10GB

library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(cowplot)

large.intestine <- readRDS("/path/to/seuratObj.rds")


###find markers for each cluster
tt.markers <- FindAllMarkers(object = large.intestine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 top50 <- tt.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
#pdf("MarkersHeatmap.pdf", width=10, height=8)
DoHeatmap(object = large.intestine, features = top10$gene, size=3) + NoLegend()
#dev.off()

#DCLK1+_progenitor_cell:
DCLK1p <- c(
"DCLK1",
"LUM",
"DCN",
"FBLN1",
"COL3A1",
"C1S",
"COL1A2",
"CXCL1",
"SPARC"
)
DCLK1p <- unique(DCLK1p)
FeaturePlot(large.intestine, features = "CD79A" , min.cutoff = "q9")
DotPlot(large.intestine, features = rev(DCLK1p), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

#LGR5+_stem_cell:
LGR5p <- c(
"LGR5",
"GDF15",
"SFN",
"PPP1R1B",
"LEFTY1",
"CHP2",
"PKP2",
"TSPAN8",
"F2RL1",
"GPX2",
"SLPI",
"SOX9",
"NOX1"
)
LGR5p <- unique(LGR5p)
DotPlot(large.intestine, features = rev(LGR5p), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

#Enterocyte:
Enter<-c(
"CDHR5",
"PRSS3",
"CYP3A5",
"EPS8L3",
"AOC1",
"DGAT1",
"NBEAL1",
"SPIB",
"MALL",
"AQP8",
"GUCA2A",
"GUCA2B",
"PRAP1",
"CEACAM7",
"LYPD8",
"PRSS3",
"CDHR5",
"MALL",
"AOC1"
)
Enter<-unique(Enter)
DotPlot(large.intestine, features = rev(Enter), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

Goblet_cell<-c(
"IFI27",
"TSPAN1",
"HLA-F",
"ERN2",
"ATP2C2",
"MUC5B",
"TSPAN1",
"CEACAM5",
"GPA33",
"CEACAM6",
"ENTPD8",
"DHRS9",
"MLPH",
"BCAS1",
"SPDEF",
"SCNN1A",
"ZG16",
"MUC5B",
"MUC4",
"ENTPD8"
)
Goblet_cell <- unique(Goblet_cell)
DotPlot(large.intestine, features = rev(Goblet_cell), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

Gobletpc<-c(
"ARL14",
"FOXA3",
"ARL14",
"SPINK4",
"WFDC2",
"REG4",
"KIAA1324",
"FOXA3",
"ZG16B",
"GALNT5",
"SLC39A7",
"TSPAN15"
)
Gobletpc <- unique(Gobletpc)
DotPlot(large.intestine, features = rev(Gobletpc), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

MKI67<-c(
"PLVAP",
"CLDN5",
"GNG11",
"RAMP2",
"TM4SF1",
"VWF",
"EGFL7",
"RAMP3",
"ENG",
"JAM2",
"ESAM",
"ECSCR",
"CLEC14A",
"RGCC",
"RHOH",
"HIST1H4C",
"GIMAP7",
"RRM2",
"BIRC5",
"UBE2C"
)
MKI67 <- unique(MKI67)
DotPlot(large.intestine, features = rev(MKI67), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

path<-c(
"SRGN",
"FYB",
"CXCR4",
"ARHGDIB",
"PTPRC",
"MAF",
"CD74",
"CCL4",
"HCST",
"DUSP2",
"RUNX3",
"PTPRC",
"GPR65",
"CD83",
"CD74",
"HLA-DRB1",
"HLA-DQA1",
"HLA-DQB1",
"HLA-DPB1",
"CXCR4",
"CD37",
"LAPTM5",
"IL1B",
"FCER1G",
"AIF1",
"HCAR3",
"MS4A6A",
"LST1",
"IL1RN",
"CLEC10A",
"MS4A7",
"CPVL",
"FAM26F",
"CLEC7A",
"CSF1R",
"CFP",
"CD209",
"CD163",
"CSF3R",
"NCF2",
"HCK",
"HCAR2",
"TLR2"
)
path <- unique(path)
DotPlot(large.intestine, features = rev(path), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

tuft<-c(
"GNG13",
"NREP",
"RGS2",
"POU2F3",
"RAC2",
"PTGS1",
"IRF7",
"FFAR3",
"ALOX5",
"TSLP",
"CD14",
"EPCAM",
"DCLK1",
"PTPRC"
)
tuft <- unique(tuft)
DotPlot(large.intestine, features = rev(tuft), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()



tcells <- c("IL7R","CCR7","TCF7","GZMK","CD8A","NKG7","CD14","FCGR3A")
DotPlot(large.intestine, features = rev(tcells), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

mono.markers<- c("CD14","CCL2","S100A9","FCGR3A","VMO1","CD86","HLA-DRA","HLA-DRB1","ITGAX")
DotPlot(large.intestine, features = rev(mono.markers), cols=c("darkblue","red"), dot.scale = 8) + RotatedAxis()

Idents(idName)="seurat_clusters"
idName <- RenameIdents(idName,
`0` = "plasma-1",
`1` = "Memory CD4+ T",
`2` = "Naive B",
`3` = "CD8+ T",
`4` = "NK/CD8+ T",
`5` = "Naive CD4+ T",
`6` = "Enterocyte progenitor",
`7` = "Goblet progenitor",
`8` = "Naive B",
`9` = "Dead",
`10` = "plasma-2",
`11` = "Paneth",
`12` = "DCLK1+ progenitor",
`13` = "Enterocyte",
`14` = "plasma-3",
`15` = "Monocytes/DC",
`16` = "MKI67+ progenitor",
`17` = "DCLK1+ progenitor",
`18` = "Paneth",
`19` = "Goblet",
`20` = "LGR5+ stem",
`21` = "mast cell",
`22` = "tuft",
`23` = "DCLK1+ tuft",
`24` = "Memory B"
)
idName$celltype.raw <- Idents(idName)

Idents(idName)="seurat_clusters"
idName <- RenameIdents(idName,
`0` = "PLA1",
`1` = "MCD4T",
`2` = "NaiB",
`3` = "CD8T",
`4` = "NKT",
`5` = "NCD4T",
`6` = "EP",
`7` = "GP",
`8` = "NaiB",
`9` = "Dead",
`10` = "PLA2",
`11` = "PAN",
`12` = "DCLKP",
`13` = "ENT",
`14` = "PLA3",
`15` = "MDC",
`16` = "MKIP",
`17` = "DCLKP",
`18` = "PAN",
`19` = "GOB",
`20` = "LGRS",
`21` = "MC",
`22` = "TUFT",
`23` = "DCLKT",
`24` = "MemB"
)
idName$celltype.short <- Idents(idName)



######  Condition markers -- DEGs ###
for(index1 in levels(idName$celltype.short)){
  index1
  subD <- subset(idName, idents = index1)
  a.markers <- FindMarkers(object = subD, ident.1 = "disease", ident.2="healthy.ctrl", group.by = "stim", verbose = FALSE)
  b.markers <- FindMarkers(object = subD, ident.1 = "disease", ident.2="self.ctrl", group.by = "stim", verbose = FALSE)
  c.markers <- FindMarkers(object = subD, ident.1 = "self.ctrl", ident.2="healthy.ctrl", group.by = "stim", verbose = FALSE)
  name1 <- paste("./DEGs/UCvHC.",index1,".csv",sep="")
  name2 <- paste("./DEGs/UCvSC.",index1,".csv",sep="")
  name3 <- paste("./DEGs/SCvHC.",index1,".csv",sep="")
  write.csv(a.markers,file=name1, quote=F)
  write.csv(b.markers,file=name2, quote=F)
  write.csv(c.markers,file=name3, quote=F)
 }
#####
