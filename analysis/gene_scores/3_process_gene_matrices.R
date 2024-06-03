library(Seurat)
library(Signac)
library(SnapATAC)
library(ArchR)
library(data.table)
library(biomaRt)
library(DelayedArray)

out_path <- "/home/siluo/public/SiyuanLuo/projects/rebuttal/gene_scores/outputs"

dataset <- "/PBMC_multiomics/"
rna <- readRDS("/home/siluo/public/SiyuanLuo/projects/sc_chromatin_1_first_results_with_2_datasets/PBMC_multi-omics_10X/rna/pbmc_rna3_with_curated_annotation.rds")

# dataset <- "/Chen_2019/"
# rna <- readRDS("/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/Chen_2019/snare_rna_3.rds")


DefaultAssay(rna) <- "RNA"
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- NormalizeData(rna, normalization.method = "LogNormalize")
rna <- ScaleData(rna)
gene_exp <- GetAssayData(object = rna, assay = "RNA", slot = "data")
hvg <- VariableFeatures(rna, assay = "RNA")


signac <- readRDS(paste0(out_path, dataset, "signac.RDS"))
archr <- readRDS(paste0(out_path, dataset, "archr.RDS"))
snapatac <- readRDS(paste0(out_path, dataset, "snapatac.RDS"))
snapatac2 <- data.frame(fread(paste0(out_path, dataset, "snapatac2.tsv")), row.names=1)
snapatac <- t(snapatac)
snapatac2 <- t(snapatac2)

genes <- intersect(intersect(intersect(intersect(hvg, rownames(signac)), rownames(archr)), rownames(snapatac)), rownames(snapatac2))


nfeatures <- 3000
feature_variance <- rna@assays$RNA@meta.features
feature_variance <- feature_variance[genes,]
feature_variance <- feature_variance[order(feature_variance$vst.variance.standardized, decreasing = TRUE)[1:1000],]

genes <- rownames(feature_variance)

newnames <- sapply(colnames(archr), function(x){strsplit(x, "#", fixed=TRUE)[[1]][2]})
names(newnames) <- NULL
colnames(archr) <- newnames

cells <- intersect(intersect(intersect(intersect(colnames(snapatac2), colnames(signac)), colnames(archr)), colnames(snapatac)), colnames(gene_exp))

# make sure each gene activity matrices have the same rows and columns
signac <- signac[genes, cells]
archr <- archr[genes, cells]
snapatac <- snapatac[genes, cells]
snapatac2 <- snapatac2[genes, cells]
gene_exp <- gene_exp[genes, cells]

saveRDS(list(signac=signac, archr=archr, snapatac=snapatac, snapatac2=snapatac2, gene_exp=gene_exp), file=paste0(out_path, dataset, "all_gene_scores.RDS"))

####################################
# get meta cell matrices

metacell_ls <- readRDS(paste0(out_path, dataset, "meta-cells.RDS"))
gene_score_ls <- readRDS(paste0(out_path, dataset, "all_gene_scores.RDS"))

gene_score_names <- c("signac", "archr", "snapatac", "snapatac2", "gene_exp")
new_list <- list()
k <- metacell_ls$idx
for(i in 1:5){
    mname <- gene_score_names[i]
    x <- as.matrix(gene_score_ls[[mname]])
    cells <- intersect(colnames(x), metacell_ls$cells)
    x <- x[, cells]
    k <- k[cells]
    x <- colsum(x, k)
    new_list[[mname]] <- x
}
saveRDS(new_list, file=paste0(out_path, dataset, "metacell_gene_scores.RDS"))