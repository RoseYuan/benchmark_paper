
import snapatac2 as snap
import scanpy as sc
import pandas as pd
import argparse


def gene_score_snapatac2(obj_file, gene_mx_file, genome):
    data = snap.read(obj_file)
    chrom_size_database = {"hg19": snap.genome.hg19,
                           "hg38": snap.genome.hg38,
                           "mm10": snap.genome.mm10}
    chrom_size = chrom_size_database[genome]
    gene_matrix = snap.pp.make_gene_matrix(data, chrom_size)

    # use the MAGIC algorithm to perform imputation and data smoothing
    sc.pp.filter_genes(gene_matrix, min_cells= 5)
    sc.pp.normalize_total(gene_matrix)
    sc.pp.log1p(gene_matrix)
    sc.external.pp.magic(gene_matrix, solver="approximate")

    df = pd.DataFrame(gene_matrix.X, index=data.obs_names, columns=gene_matrix.var_names)
    data.close()
    
    df.to_csv(gene_mx_file, sep='\t',header=True)
    return

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_obj_file", help="Full path to the input anndata file.", type=str)
    parser.add_argument("-o", "--output_gene_mx_file", help="Full path to the output file storing gene score matrix.", type=str)
    parser.add_argument("-g", "--genome", help="Genome name.", type=str)
    args = parser.parse_args()

    obj_file = args.input_obj_file
    gene_mx_file = args.output_gene_mx_file
    genome = args.genome
    
    gene_score_snapatac2(obj_file, gene_mx_file, genome)

file = "/home/siluo/public/SiyuanLuo/projects/rebuttal/number_of_features/outputs/nfeatures25k/PBMC_multiomics/feature_engineering/python/SnapATAC2/default/500/cosineCellinFile1.h5ad"

