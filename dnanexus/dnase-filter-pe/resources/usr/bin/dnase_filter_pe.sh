#!/bin/bash -e

if [ $# -ne 3 ]; then
    echo "usage v1: dnase_filter_pe.sh <unfiltered.bam> <map_threshold> <ncpus>"
    echo "Filters paired-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires samtools on path."
    exit -1; 
fi
unfiltered_bam=$1  # unfiltered bam file.
map_thresh=$2      # Mapping threshold (e.g. 10)
ncpus=$3           # Number of cpus available

unfiltered_bam_root=${unfiltered_bam%.bam}
filtered_bam_root="${unfiltered_bam_root}_filtered"
echo "-- Filtered alignments file will be: '${filtered_bam_root}.bam'"

echo "-- Filter on threashold..."
# -F 1804 means not: 0111 0000 1100
#       4 read unmapped
#       8 mate unmapped
#     256 not primary alignment
#     512 read fails platform/vendor quality checks
#    1024 read is PCR or optical duplicate
# -F 780 means:  0011 0000 1100 not: 4,8,256,512
set -x
samtools view -F 780 -f 2 -q ${map_thresh} -u $unfiltered_bam | \
        samtools sort -@ $ncpus -m 6G -n -f - ${filtered_bam_root}_tmp.sam
samtools view -hb ${filtered_bam_root}_tmp.sam > ${filtered_bam_root}_tmp.bam
samtools fixmate -r ${filtered_bam_root}_tmp.bam -O sam - | \
        samtools view -F 780 -f 2 -u - | \
        samtools sort -@ $ncpus -m 6G -f - ${filtered_bam_root}.sam
samtools view -hb ${filtered_bam_root}.sam > ${filtered_bam_root}.bam
samtools index ${filtered_bam_root}.bam
rm *.sam
set +x

echo "-- Collect bam stats..."
set -x
samtools flagstat $unfiltered_bam > ${unfiltered_bam_root}_flagstat.txt
samtools flagstat ${filtered_bam_root}.bam > ${filtered_bam_root}_flagstat.txt
samtools stats ${filtered_bam_root}.bam > ${filtered_bam_root}_samstats.txt
head -3 ${filtered_bam_root}_samstats.txt
grep ^SN ${filtered_bam_root}_samstats.txt | cut -f 2- > ${filtered_bam_root}_samstats_summary.txt
set +x

echo "-- The results..."
ls -l ${filtered_bam_root}* ${unfiltered_bam_root}_flagstat.txt
df -k .

