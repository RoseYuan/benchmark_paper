#!/bin/bash

out_path=/home/siluo/public/SiyuanLuo/projects/rebuttal/gene_scores/outputs
data_path=/home/siluo/public/SiyuanLuo/projects/benchmark/outputs

# dataset specific
# PBMC
dataset_path=PBMC_multiomics/PBMC_multiomics/feature_engineering
dataset=PBMC_multiomics
genome=hg38
gene_range_file=${out_path}/${dataset}/gene_range.RDS

# # Chen_2019
# dataset_path=Chen_2019/Chen_2019/feature_engineering
# dataset=Chen_2019
# genome=mm10
# gene_range_file=${out_path}/${dataset}/gene_range.RDS


# run SnapATAC2
echo "run SnapATAC2"
input_obj_file=${data_path}/${dataset_path}/python/SnapATAC2/default/500/cosineCellinFile1.h5ad
output_gene_mx_file=${out_path}/${dataset}/snapatac2.tsv
conda run --no-capture-output -n scATAC-benchmark python cal_gene_score.py -i $input_obj_file -o $output_gene_mx_file -g $genome

# # run Signac
echo "run Signac"
input_obj_file=${data_path}/${dataset_path}/R/Signac/all_cell_peaks/0/default/15.RDS
output_gene_mx_file=${out_path}/${dataset}/signac.RDS
Rscript cal_gene_score.R -m signac -i $input_obj_file -o $output_gene_mx_file -g $gene_range_file

# run ArchR
echo "run ArchR"
input_obj_file=${data_path}/${dataset_path}/R/ArchR/tiles/500/default/15.RDS
output_gene_mx_file=${out_path}/${dataset}/archr.RDS
Rscript cal_gene_score.R -m archr -i $input_obj_file -o $output_gene_mx_file -g $gene_range_file

# run SnapATAC
echo "run SnapATAC"
input_obj_file=${data_path}/${dataset_path}/R/SnapATAC1/default/5000/default/15.RDS
output_gene_mx_file=${out_path}/${dataset}/snapatac.RDS
Rscript cal_gene_score.R -m snapatac -i $input_obj_file -o $output_gene_mx_file -g $gene_range_file