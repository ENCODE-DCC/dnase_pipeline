#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: dnase_rep_corr.sh <density_a_starch> <density_b_starch> <rank_col> <corr_root>"
    echo "Compares two density starched bed files and calculates the correlation. Is independent of DX and encodeD."
    echo "Requires R,Rscript and chromCor.Rscript in the current directory."
    exit -1; 
fi
density_a_bed=$1  # Replicate A density starched bed produced by hotspot2
density_b_bed=$2  # Replicate B density starched bed produced by hotspot2
rank_col=$3       # Column in starch files that should be compared
corr_root=$4      # root name for output text file

echo "-- Running chromCor..."
set -x
Rscript chromCor.Rscript -c all -n $rank_col -s $density_a_bed $density_b_bed 2> chrmCor_stderr.txt | tee ${corr_root}.txt
cat chrmCor_stderr.txt >> ${corr_root}.txt
set +x
echo "-- ---- cat '${corr_root}.txt'..."
cat ${corr_root}.txt
echo "-- ------------------"

errs=`grep Error ${corr_root}.txt | wc -l`
if [ $errs -ge 1 ]; then
    echo "-- ERROR chromCor.Rscript: `grep Error ${corr_root}.txt`..."
    exit 1
fi

echo "-- The results..."
ls -l ${corr_root}*

