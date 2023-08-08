#!/bin/bash

# Function to select rows for the sample=string2 and "cell type"=string7 from file "input.tsv"

function select_rows() {
  local input_file="$1"
  local string6="$2"
  local string7="$3"
  string7=$(echo $string7 | tr '_' ' ')
  # Use awk to select rows that match the search strings in columns 6 and 7 (from 1), then take the column 1, remove prefix that ends with "+"
  awk -F"," -v col6="$string6" -v col7="$string7" '$6 ~ col6 && $7 ~ col7' "$input_file" | cut -d ',' -f1 | cut -d '+' -f 2-
}


# Function to extract all rows from a bed that contains cell ID in a barcode file
function extract_rows_by_barcodes() {
  # Get the arguments
  local bed_file="$1"
  local cell_ids_file="$2"
  # Use grep to extract the rows that contain cell IDs from the cell IDs file
  grep -Fwf "$cell_ids_file" "$bed_file" 

}



###################################################################
# main
# This script is used to subset tissue bed files to get bed file for selected cell type in the tissue. Input information: a meta data file contains the metadata of all desired cells.

tissue_bed_dir=$1
selected_cell_bed_dir=$2

input_file=$3  # meta table
output=$4

#################################
# test:
# select_rows "raw_data/selected_cell/new_metadata.csv" "artery_aorta_SM-C1PX3" "Vascular_Smooth_Muscle_1" > cells.tmp.txt
# extract_rows_by_barcodes "raw_data/GSM5589350_artery_aorta_SM-C1PX3_rep1_fragments_sorted_extracted.bed" "cells.tmp.txt" > "raw_data/selected_cell/"GSM5589350_artery_aorta_SM-C1PX3_rep1_fragments_sorted_extracted_Vascular_Smooth_Muscle_1.bed

# rm cells.tmp.txt

##################################
# Create an output file to store all the tissue-cluster pairs
if [ -f "$output" ] ; then
    rm "$output"
fi
# Initialize an associative array to store unique column pairs
declare -a unique_pairs

echo "++++++++++++ generate meta file for pairs ++++++++++++++++"
# Open file for reading
while IFS=',' read -r col1 col2 col3 col4 col5 col6 col7 col8; do
  col7=$(echo "$col7" | tr ' ' '_')
  
  # Combine the two columns into a pair and add it to the unique_pairs array
  pair="$col6,$col7"
  if [[ ! " ${unique_pairs[*]} " =~ " ${pair} " ]]; then
    unique_pairs+=("$pair")
    # echo $pair
    echo "$pair" >> $output
  fi
done < $input_file

# ##################################
# Loop over the unique column pairs
echo "++++++++++++ extract bed files ++++++++++++++++"
while IFS=',' read -r col1 col2; do
  echo $col1, $col2
  select_rows $input_file $col1 $col2 > ${selected_cell_bed_dir}/${col1}_${col2}_cells.id.txt

  file=$(find "${tissue_bed_dir}/" -name "*${col1}*_fragments_sorted_extracted.5c.bed")
  
  if [ ! -z "$file" ] # check if it's empty
  then
     extract_rows_by_barcodes $file ${selected_cell_bed_dir}/${col1}_${col2}_cells.id.txt > "${selected_cell_bed_dir}/"${col1}_fragments_sorted_extracted_${col2}.bed

     # prefix="${col1}:${col2}+"
     awk -v prefix="$col1+$col2+" 'BEGIN { FS="\t"; OFS="\t"; ORS="\r\n" } {$4=prefix$4; print}' "${selected_cell_bed_dir}/"${col1}_fragments_sorted_extracted_${col2}.bed > "${selected_cell_bed_dir}/"${col1}_fragments_sorted_extracted_${col2}_idx.bed
    # awk -v prefix="$col1:$col2+" '{$4=prefix$4; print}'
    
     cat "${selected_cell_bed_dir}/"${col1}_fragments_sorted_extracted_${col2}_idx.bed >> "${selected_cell_bed_dir}/"pooled_fragments.bed
     # # bgzip
     # bgzip -k  "${selected_cell_bed_dir}/"${col1}_fragments_sorted_extracted_${col2}.bed
     # # tabix index
     # tabix -p bed "${selected_cell_bed_dir}/"${col1}_fragments_sorted_extracted_${col2}.bed.gz
  fi
  # rm -f cell.tmp.txt
done < $output

echo "++++++++++++ sorting, compressing and indexing ++++++++++++++++"
# sort according to coordinate
sort -k1,1 -k2,2n -k3,3n "${selected_cell_bed_dir}/"pooled_fragments.bed > "${selected_cell_bed_dir}/"pooled_fragments_sorted.bed
# bgzip
bgzip "${selected_cell_bed_dir}/"pooled_fragments_sorted.bed
# tabix index
tabix -p bed "${selected_cell_bed_dir}/"pooled_fragments_sorted.bed.gz