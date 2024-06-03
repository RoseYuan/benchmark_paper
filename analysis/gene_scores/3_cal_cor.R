library(Seurat)
library(Signac)
library(SnapATAC)
library(ArchR)
library(data.table)
library(biomaRt)

out_path <- "/home/siluo/public/SiyuanLuo/projects/rebuttal/gene_scores/outputs"

# rna <- readRDS("/home/siluo/public/SiyuanLuo/projects/sc_chromatin_1_first_results_with_2_datasets/PBMC_multi-omics_10X/rna/pbmc_rna3_with_curated_annotation.rds")

dataset <- "/PBMC_multiomics/"

# dataset <- "/Chen_2019/"

# rna <- readRDS("/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/Chen_2019/snare_rna_3.rds")

# # gene_exp <- GetAssayData(object = rna, assay = "RNA", slot = "data")
# # # select highly variable genes according to the RNA-seq data
# # hvg <- VariableFeatures(rna, assay = "RNA")

# # gene_exp <- GetAssayData(object = rna, assay = "SCT", slot = "data")
# # hvg <- VariableFeatures(rna, assay = "SCT")

# # nfeatures <- 3000
# # feature.attr <- SCTResults(object = rna, slot = "feature.attributes")
# # feature_variance <- feature.attr[order(feature.attr$residual_variance, decreasing = TRUE)[1:nfeatures],]


# DefaultAssay(rna) <- "RNA"
# rna <- FindVariableFeatures(rna, nfeatures = 3000)
# rna <- NormalizeData(rna, normalization.method = "LogNormalize")
# rna <- ScaleData(rna)
# gene_exp <- GetAssayData(object = rna, assay = "RNA", slot = "data")
# hvg <- VariableFeatures(rna, assay = "RNA")


# signac <- readRDS(paste0(out_path, dataset, "signac.RDS"))
# archr <- readRDS(paste0(out_path, dataset, "archr.RDS"))
# snapatac <- readRDS(paste0(out_path, dataset, "snapatac.RDS"))
# snapatac2 <- data.frame(fread(paste0(out_path, dataset, "snapatac2.tsv")), row.names=1)
# snapatac <- t(snapatac)
# snapatac2 <- t(snapatac2)

# genes <- intersect(intersect(intersect(intersect(hvg, rownames(signac)), rownames(archr)), rownames(snapatac)), rownames(snapatac2))

# # feature_variance <- feature_variance[genes,]
# # feature_variance <- feature_variance[order(feature_variance$residual_variance, decreasing = TRUE)[1:1000],]
# nfeatures <- 3000
# feature_variance <- rna@assays$RNA@meta.features
# feature_variance <- feature_variance[genes,]
# feature_variance <- feature_variance[order(feature_variance$vst.variance.standardized, decreasing = TRUE)[1:1000],]

# genes <- rownames(feature_variance)

# newnames <- sapply(colnames(archr), function(x){strsplit(x, "#", fixed=TRUE)[[1]][2]})
# names(newnames) <- NULL
# colnames(archr) <- newnames

# cells <- intersect(intersect(intersect(intersect(colnames(snapatac2), colnames(signac)), colnames(archr)), colnames(snapatac)), colnames(gene_exp))

# # make sure each gene activity matrices have the same rows and columns
# signac <- signac[genes, cells]
# archr <- archr[genes, cells]
# snapatac <- snapatac[genes, cells]
# snapatac2 <- snapatac2[genes, cells]
# gene_exp <- gene_exp[genes, cells]

# saveRDS(list(signac=signac, archr=archr, snapatac=snapatac, snapatac2=snapatac2, gene_exp=gene_exp), file=paste0(out_path, dataset, "all_gene_scores.RDS"))


################################################################################

ls <- readRDS(file=paste0(out_path, dataset, "all_gene_scores.RDS"))
# ls <- readRDS(file=paste0(out_path, dataset, "metacell_gene_scores.RDS"))

signac <- ls$signac
archr <- ls$archr
snapatac <- ls$snapatac
snapatac2 <- ls$snapatac2
gene_exp <- ls$gene_exp

method_names <- c("Signac", "ArchR", "SnapATAC", "SnapATAC2")
gene_score_ls <- c("signac", "archr", "snapatac", "snapatac2")
output_path <- out_path


cal_cor <- function(matrix1, matrix2, method){
    # Ensure the matrices have the same dimensions
    if (all(dim(matrix1) == dim(matrix2))) {
    
    # Initialize a vector to store the correlations
    correlations <- numeric(ncol(matrix1))
    
    # Loop through the columns
    for (i in seq_len(ncol(matrix1))) {
        correlations[i] <- cor(matrix1[, i], matrix2[, i], method=method)
    }
    
    # Print the correlations
    return(correlations)

    } else {
    stop("Matrices do not have the same dimensions.")
    }
}


message("Correlation across cells:")

# for(cors in c("pearson", "spearman", "kendall")){
for(cors in c("kendall")){
    df_cor_cells <- data.frame(method=c(), value=c(), cor=c())
    for(i in 1:4){
        message(method_names[i])
        print(max(get(gene_score_ls[i])))
        print(max(gene_exp))
        value <- cal_cor(as.matrix(gene_exp), as.matrix(get(gene_score_ls[i])), method = cors)
        df_cor_cells <- rbind(df_cor_cells, data.frame(method=rep(method_names[i], length(value)), value=value, cor=rep(cors, length(value))))
    }
    # saveRDS(df_cor_cells, file = paste0(output_path, dataset, "metacells/cells_cor_df_", cors, ".RDS"))
    saveRDS(df_cor_cells, file = paste0(output_path, dataset, "cells_cor_df_", cors, ".RDS"))
}


message("Correlation across genes:")
# for(cors in c("pearson", "spearman", "kendall")){
for(cors in c("kendall")){
    df_cor_genes <- data.frame(method=c(), value=c(), cor=c())
    for(i in 1:4){
        message(method_names[i])
        print(max(get(gene_score_ls[i])))
        print(max(gene_exp))
        value <- cal_cor(as.matrix(t(gene_exp)), as.matrix(t(get(gene_score_ls[i]))), method = cors)
        df_cor_genes <- rbind(df_cor_genes, data.frame(method=rep(method_names[i], length(value)), value=value, cor=rep(cors, length(value))))
    }
    # saveRDS(df_cor_genes, file = paste0(output_path, dataset, "metacells/genes_cor_df_", cors, ".RDS"))
    saveRDS(df_cor_genes, file = paste0(output_path, dataset, "genes_cor_df_", cors, ".RDS"))
}

