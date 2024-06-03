#!/bin/bash

cd /home/siluo/public/SiyuanLuo/projects/rebuttal/integration/scripts
# conda activate integration

methods="ArchR_peaks ArchR_tiles Signac_all_cell_peaks Signac_by_cluster_peaks SnapATAC SnapATAC2_cosine SnapATAC2_jaccard aggregation"



ms=($methods)
count=${#ms[@]}
# dataset="10XPBMC"
dataset="Chen2019"

echo ""
echo ">>>>>>>>>>> Load data from previous files >>>>>>>>>>>"
for j in $(seq 1 $count); do
  m=${ms[$j-1]}
  echo $m
  python load_from_files.py -c /home/siluo/public/SiyuanLuo/projects/rebuttal/integration/scripts/configs/${dataset}_${m}.ini
done

echo ""
echo ">>>>>>>>>>> Prepare data for GLUE >>>>>>>>>>>"
for j in $(seq 1 $count); do
  m=${ms[$j-1]}
  echo $m
  python run_glue_prepare.py -c /home/siluo/public/SiyuanLuo/projects/rebuttal/integration/scripts/configs/${dataset}_${m}.ini
done

echo ""
echo ">>>>>>>>>>> Run GLUE >>>>>>>>>>>"
for j in $(seq 1 $count); do
  m=${ms[$j-1]}
  echo $m
  python run_glue.py -c /home/siluo/public/SiyuanLuo/projects/rebuttal/integration/scripts/configs/${dataset}_${m}.ini
done

echo ""
echo ">>>>>>>>>>> Evaluate GLUE >>>>>>>>>>>"
for j in $(seq 1 $count); do
  m=${ms[$j-1]}
  echo $m
  python run_glue_evaluation.py -c /home/siluo/public/SiyuanLuo/projects/rebuttal/integration/scripts/configs/${dataset}_${m}.ini
  # kernprof -l -v run_glue_evaluation.py -c /home/siluo/public/SiyuanLuo/projects/rebuttal/integration/scripts/configs/${dataset}_${m}.ini
done

# bash main_glue.sh 2>&1 | tee ./logs/2024.4.18.log