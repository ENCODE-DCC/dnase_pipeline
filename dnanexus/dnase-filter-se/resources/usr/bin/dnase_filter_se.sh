#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: dnase_filter_se.sh <unfiltered.bam> <map_threshold> <ncpus> <umi>"
    echo "Filters single-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires samtools on path."
    exit -1; 
fi
unfiltered_bam=$1  # unfiltered bam file.
map_thresh=$2      # Mapping threshold (e.g. 10)
ncpus=$3           # Number of cpus available
umi=$4             # Whether reads in bam contain UMI ids (only 'yes' means yes).

unfiltered_bam_root=${unfiltered_bam%.bam}
filtered_bam_root="${unfiltered_bam_root}_filtered"
echo "-- Filtered alignments file will be: '${filtered_bam_root}.bam'"

#    1 read paired
#    2 read mapped in proper pair
#    4 read unmapped
#    8 mate unmapped
#   16 read reverse strand
#   32 mate reverse strand
#   64 first in pair
#  128 second in pair
#  256 not primary alignment
#  512 read fails platform/vendor quality checks
# 1024 read is PCR or optical duplicate
# 2048 supplementary alignment

if [ "$umi" == "yes" ] || [ "$umi" == "y" ] || [ "$umi" == "true" ] || [ "$umi" == "t" ] || [ "$umi" == "umi" ]; then
    echo "-- UMI filtering will be performed."
    flagged_root="${unfiltered_bam_root}_post_umi"
    flagged_file="${flagged_root}.sam"
    filter_flags=`expr 512 + 8 + 4`
    # Don't think sorting by name is needed for se
    #echo "-- Sort bam by name."
    #set -x
    #samtools sort -@ $ncpus -m 6G -n -O sam -T sorted $unfiltered_bam > sorted.sam
    #set +x
    echo "-- Detect UMI with filter_reads.py."
    # NOTE script written for python3 works just as well for python2.7 as long as pysam works
    set -x
    python2.7 /usr/bin/filter_reads.py $unfiltered_bam ${flagged_root}.sam
    set +x
else
    echo "-- non-UMI filtering will be performed."
    flagged_file=$unfiltered_bam
    filter_flags=`expr 1024 + 512 + 8 + 4`
fi
 
echo "-- Filter on flags and threashold..."
# Simple version from Richard
#samtools view -F $filter_flags -b $flagged_file > ${filtered_bam_root}.bam

# More complex version
set -x
samtools view -F $filter_flags -q ${map_thresh} -u $flagged_file | \
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

