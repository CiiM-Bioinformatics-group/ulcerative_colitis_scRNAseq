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

ReadBowen <- function(bdir,bid,btype){
	bdata <- Read10X(data.dir = bdir)
	bsample <- CreateSeuratObject(counts = bdata, min.cells = 5, min.features = 200, project = "10X")
	#bsample$individual <- individual
	bsample$sampid <- bid
	bsample$samptype <- btype
	return(bsample)
}

######## UT ########
id1t <- ReadBowen(bdir="./UC_1/filtered_feature_bc_matrix",bid="id1T",btype="UC")
id3t <- ReadBowen(bdir="./UC_3/filtered_feature_bc_matrix",bid="id3T",btype="UC")
id4t <- ReadBowen(bdir="./UC_4/filtered_feature_bc_matrix",bid="id4T",btype="UC")
id5t <- ReadBowen(bdir="./UC_5/filtered_feature_bc_matrix",bid="id5T",btype="UC")
groupT4 = merge(x=id1t, y=c(id3t,id4t,id5t),project = "10X", add.cell.ids = c("id1t","id3t","id4t","id5t"))
#groupT3 = merge(x=id3t, y=c(id4t,id5t),project = "10X", add.cell.ids = c("id3t","id4t","id5t"))
######## UC ########
id2C <- ReadBowen(bdir="./SC_2/filtered_feature_bc_matrix",bid="id2C",btype="UCC")
id3C <- ReadBowen(bdir="./SC_3/filtered_feature_bc_matrix",bid="id3C",btype="UCC")
id4C <- ReadBowen(bdir="./SC_4/filtered_feature_bc_matrix",bid="id4C",btype="UCC")
id5C <- ReadBowen(bdir="./SC_5/filtered_feature_bc_matrix",bid="id5C",btype="UCC")
groupUC4 = merge(x=id2C, y=c(id3C,id4C,id5C),project = "10X", add.cell.ids = c("id2C","id3C","id4C","id5C"))
#groupUC3 = merge(x=id3C, y=c(id4C,id5C),project = "10X", add.cell.ids = c("id3C","id4C","id5C"))
######## HC ########
id6C <- ReadBowen(bdir="./HC_6/filtered_feature_bc_matrix",bid="id6C",btype="HC")
id7C <- ReadBowen(bdir="./HC_7/filtered_feature_bc_matrix",bid="id7C",btype="HC")
id8C <- ReadBowen(bdir="./HC_8/filtered_feature_bc_matrix",bid="id8C",btype="HC")
id9C <- ReadBowen(bdir="./HC_9/filtered_feature_bc_matrix",bid="id9C",btype="HC")
groupHC4 = merge(x=id6C, y=c(id7C,id8C,id9C),project = "10X", add.cell.ids = c("id6C","id7C","id8C","id9C"))


rm(id1t,
id3t,
id4t,
id5t)
rm(id2C,
id3C,
id4C,
id5C,
id6C,
id7C,
id8C,
id9C)


groupT4$stim <- "disease"
groupT4[["percent.mt"]] <- PercentageFeatureSet(object = groupT4, pattern = "^MT-")
groupT4 <- subset(x = groupT4, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25)
groupT4<-SCTransform(groupT4, vars.to.regress = "percent.mt", verbose = FALSE)


groupHC4$stim <- "healthy.ctrl"
groupHC4[["percent.mt"]] <- PercentageFeatureSet(object = groupHC4, pattern = "^MT-")
groupHC4 <- subset(x = groupHC4, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25)
groupHC4<-SCTransform(groupHC4, vars.to.regress = "percent.mt", verbose = FALSE)


groupUC4$stim <- "self.ctrl"
groupUC4[["percent.mt"]] <- PercentageFeatureSet(object = groupUC4, pattern = "^MT-")
groupUC4 <- subset(x = groupUC4, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25)
groupUC4<-SCTransform(groupUC4, vars.to.regress = "percent.mt", verbose = FALSE)

#save.image("SCgroupAll.RData")

panlist <- list(groupT4, groupHC4, groupUC4)

uc.features <- SelectIntegrationFeatures(object.list = panlist, nfeatures = 2000)
panlist <- PrepSCTIntegration(object.list = panlist, anchor.features = uc.features, verbose = FALSE)

uc.anchors <- FindIntegrationAnchors(object.list = panlist, normalization.method = "SCT", 
    anchor.features = uc.features, verbose = FALSE)
uc.intergrated <- IntegrateData(anchorset = uc.anchors, normalization.method = "SCT", 
    verbose = FALSE)
uc.intergrated <- RunPCA(uc.intergrated, verbose = FALSE)
uc.intergrated <- RunUMAP(uc.intergrated, dims = 1:30)
uc.intergrated <- FindNeighbors(object = uc.intergrated, reduction = "pca", dims = 1:30) #used to be 18
uc.intergrated <- FindClusters(uc.intergrated, resolution = 0.9)
saveRDS(uc.intergrated,file="SCgroupAll.cluster.rds")


DimPlot(object = uc.intergrated, reduction = "umap", label = TRUE) 

#### try different parameters ####
test<-uc.intergrated
DefaultAssay(test)<-"integrated"
test <- RunUMAP(test, dims = 1:18)
test <- FindNeighbors(object = test, reduction = "pca", dims = 1:18) #used to be 18
test <- FindClusters(test, resolution = 0.5)
DimPlot(object = test, reduction = "umap", label = TRUE) 
saveRDS(test,file="SCgroupAll.cluster.res05.rds")
large.intestine <- test
rm(test)

######## Standard normalization ####
DefaultAssay(large.intestine)<-"RNA"

large.intestine <- NormalizeData(large.intestine, verbose = FALSE)
large.intestine <- ScaleData(large.intestine,verbose=FALSE)


p1 <- DimPlot(object = large.intestine, reduction = "umap", group.by = "stim") 
p2 <- DimPlot(object = large.intestine, reduction = "umap", label = TRUE) 
pdf("UMAPsidebyside.pdf",width=15,height=8)
plot_grid(p1, p2)
dev.off()
pdf("UMAPbystim.pdf",width=10,height=10)
DimPlot(object = large.intestine, reduction = "umap", split.by = "sampid", ncol=3) 
dev.off()


#DimPlot(object = large.intestine, reduction = "umap", split.by = "stim", cols="red",order="disease")
p1 <- DimPlot(object = large.intestine, group.by = "stim", cols=c("grey","grey","red"),order="self.ctrl")
p2 <- DimPlot(object = large.intestine, group.by = "stim", cols=c("grey","grey","red"),order="healthy.ctrl")
p3 <- DimPlot(object = large.intestine, group.by = "stim", cols=c("grey","grey","red"),order="disease")
pdf("UMAP_3stim.pdf",width=24,height=8)
plot_grid(p1, p2, p3)
dev.off()
saveRDS(large.intestine,file="SCgroupAll.cluster.res05.rds")


