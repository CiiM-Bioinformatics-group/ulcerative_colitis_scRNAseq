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

Idents(chaoyang) = "sampid"
chaoyang <- RenameIdents(chaoyang,
`id1T`  = "Patient1_UC",
`id3T`  = "Patient3_UC",
`id4T`  = "Patient4_UC",
`id5T`  = "Patient5_UC",
`id2C`  = "Patient2_SC",
`id3C`  = "Patient3_SC",
`id4C`  = "Patient4_SC",
`id5C`  = "Patient5_SC",
`id6C`  = "Control1_HC",
`id7C`  = "Control2_HC",
`id8C`  = "Control3_HC",
`id9C`  = "Control4_HC"
)
chaoyang$sampleid <- Idents(chaoyang)

Idents(chaoyang) = "celltype.short"
chaoyang <- RenameIdents(chaoyang,
`ENT`   = "1: Enterocyte",
`EP`    = "2: Enterocyte prog",
`LGRS`  = "3: LGR5+ stem",
`GP`    = "4: Goblet prog",
`GOB`   = "5: Goblet",
`TUFT`  = "6: TRPM5+ tuft",
`DCLKT` = "7: Glial cells", ##
`DCLKP` = "8: Fibroblasts", ##
`MKIP`  = "9: Endothelial cells", ##
`PAN`   = "10: CLP/Paneth-like cells", ##
`PLA1`  = "11: Plasma-1",
`PLA2`  = "12: Plasma-2",
`PLA3`  = "13: Plasma-3",
`MemB`  = "14: Memory B",
`NaiB`  = "15: Naive B",
`NCD4T` = "16: CD4+ naive T",
`MCD4T` = "17: CD4+ memory T",
`CD8T`  = "18: CD8+ T",
`NKT`   = "19: NKT",
`MDC`   = "20: Monocytes/DC",
`MC`    = "21: Mast"
)
chaoyang$celltype.num <- Idents(chaoyang)

Idents(chaoyang) = "celltype.short"
chaoyang <- RenameIdents(chaoyang,
`ENT`   = "Enterocyte",
`EP`    = "Enterocyte prog",
`LGRS`  = "LGR5+ stem",
`GP`    = "Goblet prog",
`GOB`   = "Goblet",
`TUFT`  = "TRPM5+ tuft",
`DCLKT` = "Glial cells", ##
`DCLKP` = "Fibroblasts", ##
`MKIP`  = "Endothelial cells", ##
`PAN`   = "CLP/Paneth-like cells", ##
`PLA1`  = "Plasma-1",
`PLA2`  = "Plasma-2",
`PLA3`  = "Plasma-3",
`MemB`  = "Memory B",
`NaiB`  = "Naive B",
`NCD4T` = "CD4+ naive T",
`MCD4T` = "CD4+ memory T",
`CD8T`  = "CD8+ T",
`NKT`   = "NKT",
`MDC`   = "Monocytes/DC",
`MC`    = "Mast"
)
chaoyang$celltype.new <- Idents(chaoyang)

num.cols <- c(
"1: Enterocyte"="#92c5de", #blue
"2: Enterocyte prog"="#d1e5f0",
"3: LGR5+ stem"="#b2abd2",  #purple
"4: Goblet prog"="#fddbc7",
"5: Goblet" ="#f4a582" ,  #red
"6: TRPM5+ tuft"="#d6604d",  #
"7: Glial cells"="#8c510a",
"8: Fibroblasts"="#35978f",
"9: Endothelial cells"="#de77ae",
"10: CLP/Paneth-like cells"="#ffeda0",
"11: Plasma-1"="#0868ac",
"12: Plasma-2"="#7bccc4",
"13: Plasma-3"="#a8ddb5",
"14: Memory B"="#7fbc41",
"15: Naive B"="#276419",
"16: CD4+ naive T"="#800026",
"17: CD4+ memory T"="#cb181d",
"18: CD8+ T"="#e7298a",
"19: NKT"="#6a51a3",
"20: Monocytes/DC"="#E0B13E",
"21: Mast"="#dfc27d"
)
DimPlot(chaoyang,group.by="celltype.num",label=F,cols= num.cols)


immu.cols <- c(
"Enterocyte"="#92c5de", #blue
"Enterocyte prog"="#d1e5f0",
"LGR5+ stem"="#b2abd2",  #purple
"Goblet prog"="#fddbc7",
"Goblet" ="#f4a582" ,  #red
"TRPM5+ tuft"="#d6604d",  #
"Glial cells"="#8c510a",
"Fibroblasts"="#35978f",
"Endothelial cells"="#de77ae",
"CLP/Paneth-like cells"="#ffeda0",
"Plasma-1"="#0868ac",
"Plasma-2"="#7bccc4",
"Plasma-3"="#a8ddb5",
"Memory B"="#7fbc41",
"Naive B"="#276419",
"CD4+ naive T"="#800026",
"CD4+ memory T"="#cb181d",
"CD8+ T"="#e7298a",
"NKT"="#6a51a3",
"Monocytes/DC"="#E0B13E",
"Mast"="#dfc27d"
)


