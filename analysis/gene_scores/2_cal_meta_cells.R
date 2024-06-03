library(scDblFinder)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(RColorBrewer)

out_path <- "/home/siluo/public/SiyuanLuo/projects/rebuttal/gene_scores/outputs"
dataset <- "/PBMC_multiomics/"
# dataset <- "/Chen_2019/"

rna <- readRDS("/home/siluo/public/SiyuanLuo/projects/sc_chromatin_1_first_results_with_2_datasets/PBMC_multi-omics_10X/rna/pbmc_rna3_with_curated_annotation.rds")
# rna <- readRDS("/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/Chen_2019/snare_rna_3.rds")

DefaultAssay(rna) <- "RNA"
rna <- FindVariableFeatures(rna, nfeatures = 3000)
rna <- NormalizeData(rna, normalization.method = "LogNormalize")
rna <- ScaleData(rna)
gene_exp <- GetAssayData(object = rna, assay = "RNA", slot = "data")
hvg <- VariableFeatures(rna, assay = "RNA")

rna.sce <- as.SingleCellExperiment(rna)
metacells <- fastcluster(
  rna.sce,
  k = 500,
  rdname = "PCA",
  nfeatures = 1000,
  verbose = TRUE,
  returnType = c("metacells")
)

saveRDS(list(meta=metacells$meta, idx=metacells$idx, cells=Cells(rna.sce)), file=paste0(out_path, dataset, "meta-cells.RDS"))

########################
rna.sce <- runUMAP(rna.sce, n_neighbors = 15)
rna.sce <- runTSNE(rna.sce, n_neighbors = 15)

# Function to create a repeating color palette
get_repeating_palette <- function(palette_name, num_levels) {
  palette <- brewer.pal(min(num_levels, brewer.pal.info[palette_name, "maxcolors"]), palette_name)
  needed_repetitions <- ceiling(num_levels / length(palette))
  return(rep(palette, needed_repetitions)[1:num_levels])
}

# Generate a color palette for 500 levels
colors <- get_repeating_palette("Set2", 500)
rna.sce$metacell_idx <- factor(metacells$idx)
plotReducedDim(rna.sce, dimred="TSNE",
               colour_by="metacell_idx") + 
  theme(legend.position = "none") +
  scale_color_manual(values = colors)


