#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: dnase_filter_se.sh <unfiltered.bam> <map_threshold> <ncpus> <filtered_bam_root>"
    echo "Filters single-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires samtools on path."
    exit -1; 
fi
unfiltered_bam=$1  # unfiltered bam file.
map_thresh=$2      # Mapping threshold (e.g. 10)
ncpus=$3           # Number of cpus available
filtered_bam_root=$4 # root name for output bam (e.g. "out" will create "out.bam" and "out_flagstat.txt") 

unfiltered_bam_root=${unfiltered_bam%.bam}
echo "-- Filtered alignments file will be: '${filtered_bam_root}.bam'"

# echo "-- Handle UMI flagging and errors with 'filter_reads.py'." # No filter-reads.py UMI handling of single-end fastqs.

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

echo "-- Running picard mark duplicates on non-UMI..."
set -x
time java -Xmx4G -jar /picard/MarkDuplicates.jar INPUT=$unfiltered_bam OUTPUT=marked.bam \
  METRICS_FILE=${filtered_bam_root}_dup_qc.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
	READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
set +x
echo "-- ------------- picard MarkDuplicates"
cat ${filtered_bam_root}_dup_qc.txt
echo "-- -------------"

# Richard Sandstrom: non-UMI flags: 512 only (again, 8 and 4 are both criteria to set 512.  we don't filter dups for non UMI reads by convention).
filter_flags=512
# TODO: Optional filter out dups on non-UMI?
#filter_flags=`expr 512 + 1024`
 
echo "-- Filter on flags and threashold..."
# Simple version from Richard
set -x
samtools view -F $filter_flags -q ${map_thresh} -b $unfiltered_bam > ${filtered_bam_root}.bam
set +x
# More complex version
#set -x
#samtools view -F $filter_flags -q ${map_thresh} -u $unfiltered_bam | \
#        samtools sort -@ $ncpus -m 6G -f - ${filtered_bam_root}.sam
#samtools view -hb ${filtered_bam_root}.sam > ${filtered_bam_root}.bam
#samtools index ${filtered_bam_root}.bam
#rm *.sam
#set +x

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

