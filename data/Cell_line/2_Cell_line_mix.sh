#!/bin/bash

# Function to extract all rows from a bed that contains cell ID in a barcode file
function extract_rows_by_barcodes() {
  # Get the arguments
  local bed_file="$1"
  local cell_ids_file="$2"
  # Use grep to extract the rows that contain cell IDs from the cell IDs file
  grep -Fwf "$cell_ids_file" "$bed_file" 

}

fragment_file=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/Cell_line_mixing/GSE162690_CellLine_HighLoading.fragments.tsv.gz
fragment_file2=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/Cell_line_mixing/GSE162690_CellLine_HighLoading.fragments.tsv
fragment_file3=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/Cell_line_mixing/GSE162690_CellLine_HighLoading.fragments.prefix.tsv
cell_ids_file=/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/Cell_line_mixing/Cell_id_high.txt
output_dir=/home/siluo/public/SiyuanLuo/projects/benchmark/cleaned_data/Cell_line_mixing

mkdir -p $output_dir
gzip -dk $fragment_file
awk -F"\t" -v OFS="\t" '{ $4 = "CellLine_HighLoading+" $4; print }' $fragment_file2 > $fragment_file3 # add prefix

# filter according to QC
extract_rows_by_barcodes $fragment_file3 $cell_ids_file > "${output_dir}"/GSE162690_CellLine_HighLoading.fragments.filtered.tsv

# sort according to coordinate
sort -k1,1 -k2,2n -k3,3n "${output_dir}"/GSE162690_CellLine_HighLoading.fragments.filtered.tsv > "${output_dir}"/GSE162690_CellLine_HighLoading.fragments.filtered.sorted.tsv


fragment_file=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/Cell_line_mixing/GSE162690_CellLine_LowLoading.fragments.tsv.gz
fragment_file2=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/Cell_line_mixing/GSE162690_CellLine_LowLoading.fragments.tsv
fragment_file3=/home/siluo/public/SiyuanLuo/projects/benchmark/raw_data/Cell_line_mixing/GSE162690_CellLine_LowLoading.fragments.prefix.tsv
cell_ids_file=/home/siluo/public/SiyuanLuo/projects/benchmark/scripts/data_cleaning/Cell_line_mixing/Cell_id_low.txt
output_dir=/home/siluo/public/SiyuanLuo/projects/benchmark/cleaned_data/Cell_line_mixing

mkdir -p $output_dir
gzip -dk $fragment_file
awk -F"\t" -v OFS="\t" '{ $4 = "CellLine_LowLoading+" $4; print }' $fragment_file2 > $fragment_file3


# filter according to QC
extract_rows_by_barcodes $fragment_file3 $cell_ids_file > "${output_dir}"/GSE162690_CellLine_LowLoading.fragments.filtered.tsv

# sort according to coordinate
sort -k1,1 -k2,2n -k3,3n "${output_dir}"/GSE162690_CellLine_LowLoading.fragments.filtered.tsv > "${output_dir}"/GSE162690_CellLine_LowLoading.fragments.filtered.sorted.tsv


# merge the two
cat "${output_dir}"/GSE162690_CellLine_HighLoading.fragments.filtered.sorted.tsv >> "${output_dir}"/GSE162690_CellLine.fragments.filtered.sorted.tsv

cat "${output_dir}"/GSE162690_CellLine_LowLoading.fragments.filtered.sorted.tsv >> "${output_dir}"/GSE162690_CellLine.fragments.filtered.sorted.tsv

sort -k1,1 -k2,2n -k3,3n "${output_dir}"/GSE162690_CellLine.fragments.filtered.sorted.tsv > "${output_dir}"/GSE162690_CellLine.fragments.filtered.sorted.sorted.tsv
# bgzip
bgzip "${output_dir}"/GSE162690_CellLine.fragments.filtered.sorted.sorted.tsv
# tabix index
tabix -p bed "${output_dir}"/GSE162690_CellLine.fragments.filtered.sorted.sorted.tsv.gz
