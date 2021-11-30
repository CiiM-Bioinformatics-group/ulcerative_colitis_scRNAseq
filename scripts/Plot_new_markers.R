library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(reshape2)
library(pryr)
library(RColorBrewer)
library(tidyr)
library(cowplot)
options(future.globals.maxSize = 10737418240) 




lineage.mk <- c("PTPRC","EPCAM","COL1A1","VWF")  # PTPRC immune cells, ECAM--epi cells, COL1A1:stromal cells, VWF: glial
FeaturePlot(object = chaoyang, features=lineage.mk, order =T)
Idents(chaoyang) <- "stim2"
exp.genes = get_expressed_genes(c("control", "disease"), chaoyang, pct = 0.10)

epi <- c("PLA2G2A", "CLCA1", "REG4", "S100A14", "ITLN1", "ELF3", "PIGR", "EPCAM", "REG1B", "REG1A", "REG3A", "FABP1", "RBP2", "SST", "FABP2", "SPINK1", "FABP6", "AGR2", "AGR3", "CLDN3", "CLDN4", "DEFA6", "DEFA5", "SPINK4", "ALDOB", "LCN2", "MUC2", "KRT8", "KRT18", "TSPAN8", "OLFM4", "GPX2", "IFI27", "PHGR1", "MT1G", "CLDN7", "KRT19", "FXYD3", "LGALS4", "FCGBP", "TFF3", "TFF1")
epi <- intersect(epi,exp.genes)
m1 <- DotPlot(chaoyang, features=epi, group.by="celltype.num",cols="RdYlBu")+RotatedAxis()+ylab("")+xlab("Epithelial markers")

#epi <- c("PLA2G2A", "CLCA1", "REG4", "S100A14", "ITLN1", "ELF3", "PIGR", "EPCAM", "KRT8", "KRT18", "TSPAN8", "IFI27", "FCGBP", "TFF3")


bcs <- c("CD79A","CD27","SDC1","CD38","CD19")
tcs <- c("SELL", "LEF1", "SC5D", "CCR7","IL7R")
cytotoxic <- c("CD8A", "CD8B", "GZMA", "KLRD1", "GZMB", "CCL5", "NKG7", "CD160")
lyps <- c(bcs,tcs,cytotoxic)
m2 <- DotPlot(chaoyang, features=lyps, group.by="celltype.num",cols="RdYlBu")+RotatedAxis()+ylab("")+xlab("Lymphoid and cytotoxicity markers")

myeloid <- c("CD14","LYZ","FCGR3A","HLA-DRA","HLA-DRB1","IL1B","C5AR1","CLEC7A","CD163","PLBD1","S100A8","S100A9","S100A12","KIT", "TPSAB1")
m3 <- DotPlot(chaoyang, features= myeloid, group.by="celltype.num",cols="RdYlBu")+RotatedAxis()+ylab("")+xlab("Myeloid markers")

sto <- c("COL1A1","COL1A2","COL6A1","COL6A2","VWF","PLVAP","CDH5","S100B")
m4 <- DotPlot(chaoyang, features= sto, group.by="celltype.num",cols="RdYlBu")+RotatedAxis()+ylab("")+xlab("Stromal markers")


tcs <- c("JUN", "KLF6", "FOSB", "PTGER2", "FOS", "SYTL3", "SPRY1", "ANKRD28", "GPR171", "PDE4D",
"JAML", "IL7R", "GLIPR1", "CD69", "NFKBIA", "PPP1R15A", "NFKBIZ", "TNFAIP3", "PTGER4", "ANXA1")

IFITM3
tcs <- intersect(tcs,exp.genes)





imm.reg <- c("TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFRSF1B", "AC133644.2", "CTLA4", "TIGIT", "ICA1", "RHBDD2", "MAGEH1", "IL2RA", "TBC1D4", "BATF", "IKZF2", "FOXP3")

DotPlot(chaoyang, features=epi, group.by="celltype.new")



HuangBcell1 <- c("CD19","CD27","CD38","IGHD","IGHM","IGJ",
"IGHG1","IGHA1","CD69","CD44","KLF2","SELL","ITGA4","CCR7","FCER2","S1PR1")
HuangBcell2 <- c("BANK1","CR2","CXCR4","CXCR5","CD40","FCRL4","ITGAE","SERPINA9",
"AICDA","PCNA","MKI67","SDC1","TNFRSF17")
DotPlot(chaoyang, features= c(HuangBcell1, HuangBcell2), group.by="celltype.num",cols="RdYlBu")+RotatedAxis()+ylab("")+xlab("B markers")



####### plot singler labels####

chaoyang$blueprint.flt <- as.character(chaoyang$blueprint.fines)
a <- table(chaoyang$blueprint.flt) < 50
a <- a[a %in% TRUE]
chaoyang$blueprint.flt[chaoyang$blueprint.flt %in% names(a)] = "other"
DimPlot(chaoyang,group.by = "blueprint.flt",label = T)

chaoyang$hpca.flt <- as.character(chaoyang$hpca.fines)
a <- table(chaoyang$hpca.flt) < 100
a <- a[a %in% TRUE]
chaoyang$hpca.flt[chaoyang$hpca.flt %in% names(a)] = "other"
DimPlot(chaoyang,group.by = "hpca.flt",label = T)



######## plot markers for immune and inst cells separately ####

inst.s <- c("GOB","GP","LGRS","ENT","EP","TUFT","DCLKT","DCLKP","MKIP","PAN")
immu.s <- c("MemB","NaiB","PLA1","PLA2","PLA3","NCD4T","MCD4T","CD8T","NKT","MDC","MC")
inst   <- c("Goblet","Goblet prog", "LGR5+ stem","Enterocyte","Enterocyte prog","TRPM5+ tuft","DCLK1+ tuft",
"DCLK1+ prog","MKI67+ prog","Paneth")
immu <- c("Memory B","Naive B","Plasma-1","Plasma-2","Plasma-3","Naive CD4+ T","Memory CD4+ T","CD8+ T","NKT","Monocytes/DC","Mast")

tt.top5 <- tt.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
tt.top5$cluster <- factor(tt.top5$cluster,levels = c(inst.s,immu.s))
tt.top5.sort <- tt.top5[order(tt.top5$cluster),]
top5.gene <- unique(tt.top5.sort$gene)


dat <- NULL
dnames <- NULL
Idents(chaoyang) <- "celltype.new"
for(index1 in c(inst,immu)){
  subD <- subset(chaoyang, idents = index1)
  sdat <- rowMeans(subD@assays$RNA@scale.data)
  dat <- cbind(dat,sdat)
  dnames <- c(dnames, index1)
 }
 colnames(dat) <- dnames
dat.top5 <- dat[top5.gene,]
pheatmap(dat.top5,scale="row",cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))

#################
chaoyang.inst <- subset(chaoyang, idents=inst)
chaoyang.immu <- subset(chaoyang, idents=immu)
#################




get.violin.data1 <- function(seurat, genes) {
  #  output = data.frame(gene = character(0), value= numeric(0), ident = character(0))
  #####################################################################################################
  output = data.frame(gene = character(0), value= numeric(0), ident = character(0), tech=character(0))
  #####################################################################################################
  for (gene in genes) {
    if(any(gene == seurat@assays$RNA@data@Dimnames[[1]])){
      data.use = data.frame(FetchData(seurat,gene))
      data.use = t(data.use)
      data.melt=data.frame(rep(gene, length(seurat@active.ident)))
      colnames(data.melt)[1]="gene"
      data.melt$value=as.numeric(data.use[1,1:length(seurat@active.ident)])
      data.melt$id=names(data.use)[1:length(seurat@active.ident)]
      data.melt$ident=seurat@active.ident
      ############################################
      #data.melt$tech=seurat$stim # ???What was used to stimulate the cells?
      ############################################
      if(any(data.melt$value != 0)) noise = rnorm(length(data.melt$value))/100000 else noise = 0
      data.melt$value=as.numeric(as.character(data.melt$value))+noise
      output = rbind(output, data.melt)
    } else {
      data.melt=data.frame(
        gene = rep(gene, seurat@assays$RNA@data@Dim[2]),
        value = rep(0, seurat@assays$RNA@data@Dim[2]),
        ident = seurat@active.ident
      )
      output = rbind(output, data.melt)
    }
  }
  return(output)
}


inst.glist <- c(
"ZG16", "SPINK4", "FCGBP",
"MUC4", "KRT19",
"GDF15", "PPP1R1B", "LEFTY1",
"CDHR5", "PRAP1","LYPD8",
"TSPAN8", "GPX2",
"PTGS1", "TRPM5",
"CRYAB","CLU","SPARC","DCLK1",
"LUM", "DCN",
"PLVAP", "CLDN5",
"RAC2","COTL1","ARHGDIB"
)
violin.plot.data <- get.violin.data1(chaoyang.inst, inst.glist)
ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) +
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=70),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")

immu.glist <- c(
"CD79A", "CD27",  "SDC1", "CD38",  "CD19",# 24
"CCR7", "TCF7","IL7R",
"GZMK","CD8A","NKG7",
"CD14", "CST3", "LYZ",
"CD69", "KIT", "TPSAB1"
)
violin.plot.data <- get.violin.data1(chaoyang.immu, immu.glist)
ggplot(violin.plot.data, aes(factor(ident),value)) +
  geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(ident)),show.legend = F) +
  ylab("") + xlab("") +
  coord_flip() +facet_wrap(~ gene,scales = "free_x", ncol = length(levels(violin.plot.data$gene))) +
  theme(strip.text.x = element_text(size=12, angle=70),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        panel.spacing.x = unit(c(-0.2), "lines"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin = unit(c(0.5,0.4,0.5,1.3), "cm")) +
  scale_x_discrete(limits = rev(levels(violin.plot.data$ident)), position = "top")


my.color <- hue_pal()(22)
grey.44 <- rep("grey",44)

chaoyang$ct.stim <- paste(chaoyang$celltype.short,chaoyang$samptype,sep=".")
aaa <- c(inst.s, immu.s, "Dead")
aaa.uc <- paste(aaa, "UC",sep=".")
aaa.hc <- paste(aaa, "HC",sep=".")
aaa.UCC <- paste(aaa, "UCC",sep=".")
DimPlot(chaoyang, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.uc)
DimPlot(chaoyang, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.hc)
DimPlot(chaoyang, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.UCC)

my.color <- hue_pal()(10)
grey.44 <- rep("grey",20)

chaoyang.inst$ct.stim <- paste(chaoyang.inst$celltype.short,chaoyang.inst$samptype,sep=".")
aaa.uc <- paste(inst.s, "UC",sep=".")
aaa.hc <- paste(inst.s, "HC",sep=".")
aaa.UCC <- paste(inst.s, "UCC",sep=".")
DimPlot(chaoyang.inst, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.uc)
DimPlot(chaoyang.inst, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.hc)
DimPlot(chaoyang.inst, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.UCC)

my.color <- hue_pal()(11)
grey.44 <- rep("grey",22)
immu.s <- c("PLA1","PLA2","PLA3","MemB","NaiB","NCD4T","MCD4T","CD8T","NKT","MDC","MC")
chaoyang.immu$ct.stim <- paste(chaoyang.immu$celltype.short, chaoyang.immu$samptype,sep=".")
aaa.uc <- paste(immu.s, "UC",sep=".")
aaa.hc <- paste(immu.s, "HC",sep=".")
aaa.UCC <- paste(immu.s, "UCC",sep=".")
DimPlot(chaoyang.immu, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.uc)  + NoLegend()
DimPlot(chaoyang.immu, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.hc)  + NoLegend()
DimPlot(chaoyang.immu, group.by="ct.stim", cols = c(grey.44,rev(my.color)), order = aaa.UCC) + NoLegend()


Idents(chaoyang) <- "samptype"
chaoyang.hc <- subset(chaoyang, idents="HC")
chaoyang.sc <- subset(chaoyang, idents="UCC")
chaoyang.uc <- subset(chaoyang, idents="UC")

dat <- NULL
dnames <- NULL
Idents(chaoyang.uc) <- "celltype.new"
Idents(chaoyang.sc) <- "celltype.new"
#for(index1 in c(inst,immu,"Dead")){
for(index1 in immu){
  subUC <- subset(chaoyang.uc, idents = index1)
  subHC <- subset(chaoyang.sc, idents = index1)
  sdat <- rowMeans(subUC@assays$RNA@data) - rowMeans(subHC@assays$RNA@data)
  dat <- cbind(dat,sdat)
  dnames <- c(dnames, index1)
 }
 colnames(dat) <- dnames

list.anti.proc <- read.table(file = "./symbol_list/list.antigen_processing.txt",stringsAsFactors = F)[,1]
list.allo.reje <- read.table(file = "./symbol_list/list.allograft.rejection.txt",stringsAsFactors = F)[,1]
list.IgA.prod <- read.table(file = "./symbol_list/list.IgA_production.txt",stringsAsFactors = F)[,1]
list.IL17.path <- read.table(file = "./symbol_list/list.IL17_pathway.txt",stringsAsFactors = F)[,1]
list.th1.path <- read.table(file = "./symbol_list/list.Th1andTh2_diff.txt",stringsAsFactors = F)[,1]

 a<- intersect(list.IL17.path,row.names(dat))
dat.il17 <- dat[a,]
new.dat.il17<- dat.il17[rowSums(abs(dat.il17))>1,]
pheatmap(new.dat.il17, breaks=seq(-3,3,0.06),cluster_rows = T,cluster_cols = T, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),treeheight_row=0)

 b<- intersect(list.allo.reje,row.names(dat))
 dat.allo <- dat[b,]
new.dat.allo <- dat.allo[rowSums(abs(dat.allo))>1,]
pheatmap(new.dat.allo, breaks=seq(-2,2,0.04),cluster_rows = T,cluster_cols = T, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),treeheight_row=0)

 a<- intersect(th17,row.names(dat))
dat.th17 <- dat[a,]
new.dat.th17 <- dat.th17[rowSums(abs(dat.th17))>0.5,]
pheatmap(new.dat.th17, breaks=seq(-2,2,0.04),cluster_rows = T,cluster_cols = T, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),treeheight_row=0)


 a<- intersect(list.IgA.prod,row.names(dat))
dat.iga<- dat[a,]
new.dat.iga <- dat.iga[rowSums(abs(dat.iga))>0.5,]
pheatmap(new.dat.iga, breaks=seq(-2,2,0.04),cluster_rows = T,cluster_cols = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),treeheight_row=0)


 a<- intersect(list.th1.path,row.names(dat))
dat.iga<- dat[a,]
new.dat.iga <- dat.iga[rowSums(abs(dat.iga))>0.5,]
pheatmap(new.dat.iga, breaks=seq(-2,2,0.04),cluster_rows = T,cluster_cols = T, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),treeheight_row=0, ,treeheight_col=20)


