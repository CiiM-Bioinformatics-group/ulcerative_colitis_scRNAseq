# ulcerative_colitis_scRNAseq
Scripts: R and Perl codes used in ulcerative_colitis_scRNAseq projects
metadata.final.txt contains the cellular information and annotations

# analysis description
scRNA-seq was performed with 10x genomics.Reads from fastq was aligned to Human genome (Hg19) by cellranger(v3.0.1).
Count matrix of each sample was stored in countTables.

Seurat v3.1 was applied to analyze the count matrix from each sample. To filter out low-quality cells and doublets, empirically filtering criteria were applied to each cell: number of estimated genes should be higher than 100 and lower than 6000 and the ratio of reads mapping to the mitochondria should be lower than 25%. Only genes detected in at least 5 cells were maintained for subsequent analyses.

Seurat v3.1 integration workflow with SCTransform normalization method was used to cluster cells from different samples into distinct cell subsets62. We followed this workflow with the following steps: Firstly, we SCTransformed each sample and merged them into UC, SC and HC datasets. Next, we selected 2,000 variable features among three datasets and identified anchors from these features to integrate the datasets. These two steps corrected batch effects and prevented cells clustering by patients or disease phenotypes rather than by cell types or cell subsets. Principal component analysis (PCA) has then been performed on the integrated datasets, followed by Shared Nearest Neighbor (SNN) Graph construction using PC1 to 20 and k=20 nearest neighbours to identify unsupervised cell clusters. Finally, Uniform Manifold Approximation and Projection (UMAP) was used to visualize the cell clusters. 

In order to keep the biological differences for downstream analyses, the above-mentioned batch correction was only used in the cell clustering and PCA related steps. For the other analyses, we used standard LogNormalization methods. The original gene counts for each cell were normalized by total UMI counts and multiplied by 10,000 (TP10K), and then log transformed by log (TP10K+1).

Seurat object was saved in Folder: "analysis"
R-scripts were saved in Folder: "scripts"

# Sample information
Totally 12 samples, including:
UC, inflamed biopsies from UC patients: id1T, id3T, id4T, id5T. 
SC, non-inflamed biopsies from UC patients: id2C, id3C, id4C, id5C.
HC, healthy biopsies from healthy individuals: id6C, id7C, id8C, id9C.

Sample IDs in meatadata.list.txt includes sample IDs and their condition, where UC stands for UC-disease-tissue, SC stands for healthy tissue from UC patient (self control), HC means healthy control (tissue from healthy patient). id1, id2, ... stands for different individual, e.g. id3C(SC) and id3T(UC) were two samples from same individual.

Sample id2T was pre-filtered due to low overall quality according to the reads mapping. Then, we performed a single-cell RNA sequencing (scRNA-seq) analysis of 12 colon biopsies from 5 UC patients including 4 inflamed (UC), 4 non-inflamed (self-control, SC) biopsies and 4 healthy biopsies (HC) from healthy individuals.

# scripts information:
loadData_and_clustering.R       # Load and clustering of the whole data (as shown in Figure 1 in the publication)
update_metadata.R  	        # update the celltype identifications and others
test_markers.R	                # test markers and DEGs (used in Figure 1)
bakup_old_celltypemarkers.txt	#celltype markers collected from cellMarkerDB http://bio-bigdata.hrbmu.edu.cn/CellMarker
foundGWAS	                #Match DE genes to GWAS reported genes, GWAS collected from gwas catalog https://www.ebi.ac.uk/gwas/
plotGWASgenes.R	                #Plot the matched GWAS reported genes (as shown in Figure 4)
Plot_new_markers.R	        #Plot the selected markers (as shown in Figure 2,3 and suppl Figures)


# Link to publication:
https://www.sciencedirect.com/science/article/pii/S2352345X21000266



