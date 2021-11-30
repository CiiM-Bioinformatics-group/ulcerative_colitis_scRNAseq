
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

gwas <- read.table("./foundGWAS/test",header=T,sep="\t")
gwas$show <- "no"
gwas$show[abs(gwas$HC)>1] <- "yes"
gwas$show[abs(gwas$SC)>1] <- "yes"
gwas$show[gwas$HC >0 & gwas$SC <0] <- "yes"

p <- ggplot(gwas, aes(x=HC, y=SC, label=Gene)) + geom_point(aes(color=Cell), size=2)
#p <- ggplot(gwas.imm, aes(x=HC, y=SC, label=Gene)) + geom_point(aes(color=factor(Cell, levels=immu.s)), size=3)

gwas.int <- gwas[gwas$Cell %in% inst.s,]
p <- ggplot(gwas.int, aes(x=HC, y=SC, label=Gene)) + geom_point(aes(color=factor(Cell, levels=inst.s)), size=3, show.legend = F)
p2 <- p +  geom_hline(aes(yintercept=0),linetype="dashed") +
    geom_vline(aes(xintercept=0),linetype="dashed") + 
    theme_bw() + xlim(c(-2.5,2.5))+ ylim(c(-2,2)) +
    geom_label_repel(aes(label = ifelse(show  == "yes", as.character(Gene), '')), hjust = 1.25, vjust = 0, size = 4, alpha=0.5,na.rm=TRUE,seed=42)+
    geom_label_repel(aes(label = ifelse(show  == "yes", as.character(Gene), '')), hjust = 1.25, vjust = 0, size = 4, alpha=1,na.rm=TRUE,seed=42,fill=NA)
    p2 + geom_abline(intercept = 0, slope=1, color="darkgrey", linetype="dashed")+ theme(text=element_text(size=20)) +
    xlab("logFoldChange (UC vs HC)") + ylab("logFoldChange (UC vs SC)")


gwas.imm <- gwas[gwas$Cell %in% immu.s,]
p <- ggplot(gwas.imm, aes(x=HC, y=SC, label=Gene)) + geom_point(aes(color=factor(Cell, levels=immu.s)), size=3, show.legend = F)
p2 <- p +  geom_hline(aes(yintercept=0),linetype="dashed") +
    geom_vline(aes(xintercept=0),linetype="dashed") + 
    theme_bw() + xlim(c(-2.5,2.5))+ ylim(c(-2,2)) +
    geom_label_repel(aes(label = ifelse(show  == "yes", as.character(Gene), '')), hjust = 1.25, vjust = 0, size = 4)
    p2 + geom_abline(intercept = 0, slope=1, color="darkgrey", linetype="dashed")+ theme(text=element_text(size=20)) +
    xlab("logFoldChange (UC vs HC)") + ylab("logFoldChange (UC vs SC)")



