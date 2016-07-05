#!/bin/bash -e

if [ $# -ne 6 ]; then
    echo "usage v1: dnase_idr.sh <peaks_a_bed> <peaks_b_bed> <format> <rank_col> <chrom_sizes> <idr_root>"
    echo "Compares two sets of peaks and generates an Irreproducible Discovery Rate report. Is independent of DX and encodeD."
    echo "If format is 'broadPeak' or 'broadPeak' then output will be converted to bigBed."
    exit -1; 
fi
peaks_a_bed=$1  # Rampage dumplicates marked bam file.
peaks_b_bed=$2  # long-RNA-seq control bam file.
format=$3       # File input format: bed, narrowPeak, broadPeak
rank_col=$4     # Column that idr should use for ranking (signal.value, q.value, score, index # for bed)
chrom_sizes=$5  # chrom_sizes file that matches the genome used to create bam_root.
idr_root=$6     # root name for output bb and png files 

peaks_a_file=$peaks_a_bed
if [[ "$peaks_a_file" == *.gz ]]; then
    echo "-- Uncompressing ${peaks_a_file} ..."
    set -x
    gunzip $peaks_a_file
    set +x
    peaks_a_file=${peaks_a_file%.gz}
fi    
peaks_b_file=$peaks_b_bed
if [[ "$peaks_b_file" == *.gz ]]; then
    echo "-- Uncompressing ${peaks_b_file} ..."
    set -x
    gunzip $peaks_b_file
    set +x
    peaks_b_file=${peaks_b_file%.gz}
fi    

echo "-- Removing any spike-ins from bed files..."
set -x
grep "^chr" $peaks_a_file | grep -v "^chrEBV" > peaks_a_clean.bed
grep "^chr" $peaks_b_file | grep -v "^chrEBV" > peaks_b_clean.bed
set -x

echo "-- Running IDR..."
set -x
idr/bin/idr --input-file-type $format --rank $rank_col --plot --verbose --samples peaks_a_clean.bed peaks_b_clean.bed 2>&1 | tee idr_summary.txt
sort -k1,1 -k2,2n < idrValues.txt > ${idr_root}.bed
mv idrValues.txt.png ${idr_root}.png
set +x

if [ "$format" == "broadPeak" ]; then 
    echo "* Converting bed to 'broad_idr' bigBed..."
    set -x
    bedToBigBed ${idr_root}.bed -type=bed6+ -as=/usr/bin/broad_idr.as $chrom_sizes ${idr_root}.bb
    set +x
elif [ "$format" == "narrowPeak" ]; then 
    echo "* Converting bed to 'narrow_idr' bigBed..."
    set -x
    bedToBigBed ${idr_root}.bed -type=bed6+ -as=/usr/bin/narrow_idr.as $chrom_sizes ${idr_root}.bb
    set +x
fi

echo "-- The results..."
ls -l ${idr_root}*

