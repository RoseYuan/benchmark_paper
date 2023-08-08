In this directory, we provide the link and code we used to download and preprocess the scATAC-seq and scRNA-seq data. At the same time, we deposited all the processed datasets on [Zenodo](https://zenodo.org/record/8212920), which can be directly input into our [Snakemake pipeline](https://github.com/RoseYuan/sc_chromatin_benchmark).

*** 
- [Downloads](#downloads)
  - [Cell line](#cell-line)
  - [Human adult atlas](#human-adult-atlas)
  - [10XPBMC](#10xpbmc)
  - [Chen2019](#chen2019)
  - [Buenrostro2018](#buenrostro2018)

*** 

# Downloads<a name="Downloads"></a>
## Cell line<a name="cellline"></a>
1. ATAC fragment files
    ```commandline
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162690/suppl/GSE162690_CellLine_HighLoading.fragments.tsv.gz
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162690/suppl/GSE162690_CellLine_LowLoading.fragments.tsv.gz
    ```
2.  Annotation files
    ```commandline
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162690/suppl/GSE162690_demuxlet_CellLineHigh_mergedBAMs.best.txt.gz
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162690/suppl/GSE162690_demuxlet_CellLineLow_mergedBAMs.best.txt.gz
    ```


## Human adult atlas<a name="atlas"></a>
Selected & downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184462.


## 10XPBMC<a name="pbmc"></a>
1. ATAC fragment files
    ```commandline
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
    ```
2. RNA
    ```commandline
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5

    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
    ```

## Chen2019<a name="chen2019"></a>
1. ATAC 
   
   The fragment file is provided by an [online tutorial](https://stuartlab.org/signac/articles/snareseq.html).
    ```commandline
    wget https://signac-objects.s3.amazonaws.com/snareseq/fragments.sort.bed.gz
    wget https://signac-objects.s3.amazonaws.com/snareseq/fragments.sort.bed.gz.tbi
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074_AdBrainCortex_SNAREseq_chromatin.barcodes.tsv.gz
    ```
2. RNA
   
   The Seurat object of the reference dataset used for label transfer is provided by an [online tutorial](https://stuartlab.org/signac/articles/snareseq).
    ```commandline
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074_AdBrainCortex_SNAREseq_cDNA.barcodes.tsv.gz
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074_AdBrainCortex_SNAREseq_cDNA.counts.mtx.gz
    wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv.gz

    # rename, necessary for function Seurat::Read10X()
    mv GSE126074_AdBrainCortex_SNAREseq_cDNA.barcodes.tsv.gz barcodes.tsv.gz
    mv GSE126074_AdBrainCortex_SNAREseq_cDNA.counts.mtx.gz matrix.mtx.gz
    mv GSE126074_AdBrainCortex_SNAREseq_cDNA.genes.tsv.gz features.tsv.gz

    # download the reference dataset
    wget https://signac-objects.s3.amazonaws.com/allen_brain.rds
    ```

## Buenrostro2018<a name="Buenrostro2018"></a>
1. ATAC bam files

    The bam files are downloaded from here: https://www.dropbox.com/sh/8o8f0xu6cvr46sm/AAB6FMIDvHqnG6h7athgcm5-a/Buenrostro_2018.tar.gz?dl=0.

    This link is provided online, see [GitHub](https://github.com/pinellolab/scATAC-benchmarking/tree/master/Real_Data/Buenrostro_2018/input/sc-bams_nodup) or [bioRxiv](https://www.biorxiv.org/content/10.1101/739011v1.full).

2. Annotation
   
    The meta file is from https://github.com/pinellolab/scATAC-benchmarking/blob/master/Real_Data/Buenrostro_2018/input/metadata.tsv.
