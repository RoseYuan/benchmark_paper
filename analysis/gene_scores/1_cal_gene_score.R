suppressPackageStartupMessages({
    library(optparse)
    library(GenomicRanges)
    })

# Signac, return a sparse matrix of class ngCMatrix
gene_score_signac <- function(obj_file, gene_mx_file){
    suppressPackageStartupMessages({
    require(Signac)
    library(SeuratObject)
    library(Seurat)
    })
    obj <- readRDS(obj_file)
    gene_score <- GeneActivity(obj,
    assay=NULL,
    extend.upstream = 2000,
    extend.downstream = 0,
    max.width = 5e+05,
    biotypes = "protein_coding")
    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    obj[['RNA']] <- CreateAssayObject(counts = gene_score)
    obj <- NormalizeData(
    object = obj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(obj$nCount_RNA)
    )
    gene_score_normalized <- GetAssayData(object = obj, assay = "RNA", slot = "data")
    saveRDS(gene_score_normalized, gene_mx_file)
}

gene_score_ArchR <- function(obj_file, gene_mx_file, gene_range_file){
    suppressPackageStartupMessages({
	require(ArchR)
    })
    obj <- readRDS(obj_file)
    # obj@sampleColData$ArrowFiles <- "/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/Chen_2019/Chen_2019/feature_engineering/R/ArchR/tiles/500/default/proj/ArrowFiles/CellinFile1.arrow"
    
    obj <- addGeneScoreMatrix(
    input = obj,
    genes = getGenes(obj),
    geneModel = "exp(-abs(x)/5000) + exp(-1)",
    matrixName = "GeneScoreMatrix",
    extendUpstream = c(1000, 1e+05),
    extendDownstream = c(1000, 1e+05),
    geneUpstream = 5000,
    geneDownstream = 0,
    useGeneBoundaries = TRUE,
    useTSS = FALSE,
    extendTSS = FALSE,
    tileSize = 500,
    ceiling = 4,
    geneScaleFactor = 5,
    scaleTo = 10000,
    force = TRUE
    )
    gene_score_obj <- getMatrixFromProject(
    ArchRProj = obj,
    useMatrix = "GeneScoreMatrix"
    )
    gene_score <- assays(gene_score_obj)$GeneScoreMatrix
    rownames(gene_score) <- rowData(gene_score_obj)$name

    saveRDS(log(gene_score+1), gene_mx_file)
    saveRDS(getGenes(obj), gene_range_file)
}

gene_score_snapatac <- function(obj_file, gene_mx_file, gene_range_file){
    suppressPackageStartupMessages({
    require(SnapATAC)
    })
    genes.gr <- readRDS(gene_range_file)
    genes.gr <- genes.gr[!is.na(genes.gr$symbol),]
    genes.gr$name <- genes.gr$symbol

    x.sp <- readRDS(obj_file)
    # x.sp@file <- rep("/home/siluo/public/SiyuanLuo/projects/benchmark/outputs/Chen_2019/Chen_2019/feature_engineering/R/SnapATAC1/default/5000/default/Chen_2019_filtered_fragments_sorted.bed.snap", 5199)
    
    # re-add the cell-by-bin matrix to the snap object;
    x.sp <- addBmatToSnap(x.sp)
    x.sp <- createGmatFromMat(
    obj=x.sp, 
    input.mat="bmat",
    genes=genes.gr,
    do.par=TRUE,
    num.cores=10)
    
    # normalize the cell-by-gene matrix
    x.sp <- scaleCountMatrix(
    obj=x.sp, 
    cov=x.sp@metaData$UQ + 1,  # the meaning of column names of metadata is here: https://github.com/r3fang/SnapTools/blob/master/README.md
    mat="gmat",
    method = "logRPM")
    
    # smooth the cell-by-gene matrix
    x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:15,
    k=15
    )
    x.sp <- runMagic(
    obj=x.sp,
    input.mat="gmat",
    step.size=3)
    
    rownames(x.sp@gmat) <- rownames(x.sp@bmat)
    saveRDS(x.sp@gmat, gene_mx_file)
}



option_list <- list(
	make_option(c("-m", "--method"), type="character", default=NA, help="method to use"),
	make_option(c("-i", "--input_obj_file"), type="character", default=NA, help="input RDS file path."),
	make_option(c("-o", "--output_gene_mx_file"), type="character", default=NA, help="output file for gene score matrix."),
	make_option(c("-g", "--gene_range_file"), type="character", default=NA, help="file storing the gene range object.")
)

opt <- parse_args(OptionParser(option_list=option_list))

method <- opt$method
obj_file <- opt$input_obj_file
gene_mx_file <- opt$output_gene_mx_file
gene_range_file <- opt$gene_range_file

if(tolower(method) == "signac") {
    gene_score_signac(obj_file, gene_mx_file)
}else if (tolower(method) == "archr") {
   gene_score_ArchR(obj_file, gene_mx_file, gene_range_file)
}else if (tolower(method) == "snapatac") {
   gene_score_snapatac(obj_file, gene_mx_file, gene_range_file)
}else{
    stop("Unrecognized method!")
}
