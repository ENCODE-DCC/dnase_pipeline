#!/bin/bash -e

if [ $# -ne 8 ]; then
    echo "usage v1: dnase_pool_reps.sh <bam_A> <bam_B> <peaks_A.bb> <peaks_B.bb> <signals_A.bw> <signals_B.bw> <chrom.sizes> <out_root>"
    echo "Pools bam, merges peaks and correlates signals from two replicates.  Is independent of DX and encodeD."
    echo "Requires samtools, bedtools, bigBedToBed, bedToBigBed, bigWigCorrelate, edwComparePeaks, edwBamStats on path."
    exit -1; 
fi
bam_A=$1        # First bam file to pool.
bam_B=$2        # Second bam file to pool.
peaks_A=$3      # First narrowPeaks bigBed file to merge.
peaks_B=$4      # Second narrowPeaks bigBed file to merge.
signals_A=$5    # First narrowPeaks bigWig file to correlate.
signals_B=$6    # Second narrowPeaks bigWig file to correlate.
chrom_sizes=$7  # Chrom sizes file that matches the reference genome which bam reads were aligned to.
out_root=$8     # root name for pooled/merge files (e.g. "out" will create "out_pooled.bam", "out_merged_narrowPeak.bb")

bam_pooled_root="${out_root}_pooled"
peaks_merged_root="${out_root}_merged_narrowPeak" 
signal_root_root="${out_root}_signal" 
echo "-- output: '${bam_pooled_root}.bam', '${bam_pooled_root}_edwBamStats.txt', "
echo "           '${peaks_merged_root}.bed/.bb', '${out_root}_peaks_overlap_qc.txt', '${signal_root}_corr_qc.txt'"

echo "-- Merging bams..."
set -x
samtools cat $bam_A $bam_A > ${bam_pooled_root}.bam
set +x

echo "-- Running edwBamStats on '${bam_pooled_root}.bam'"
set -x
edwBamStats ${bam_pooled_root}.bam ${bam_pooled_root}_edwBamStats.txt
set +x

echo "-- Merging peaks..."
set -x
bigBedToBed $peaks_A peaks_A.bed
bigBedToBed $peaks_B peaks_B.bed
cat peaks_A.bed peaks_B.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > ${peaks_merged_root}.bed
bedToBigBed ${peaks_merged_root}.bed $chrom_sizes ${peaks_merged_root}.bb
set +x

echo "-- Correlating signals..."
set -x
bigWigCorrelate -restrict=${peaks_merged_root}.bb $signals_A $signals_B > ${signal_root}_corr_qc.txt
set +x

echo "-- Comparing peak overlaps..."
set -x
edwComparePeaks $peaks_A $peaks_B ${out_root}_peaks_overlap_qc.txt
set +x

echo "-- The results..."
ls -l ${bam_pooled_root}* ${peaks_merged_root}* ${broadPeak_root}* ${signal_root}*
df -k .

