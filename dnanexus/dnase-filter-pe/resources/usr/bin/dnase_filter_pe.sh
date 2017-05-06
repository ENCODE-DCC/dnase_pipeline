#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: dnase_filter_pe.sh <unfiltered.bam> <map_threshold> <ncpus> <umi> <filtered_bam_root>"
    echo "Filters paired-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires samtools, java8 on path, and (umi aware) picard.jar and filter_reads.py in working directory."
    exit -1; 
fi
unfiltered_bam=$1  # unfiltered bam file.
map_thresh=$2      # Mapping threshold (e.g. 10)
ncpus=$3           # Number of cpus available
umi=$4             # Whether reads in bam contain UMI ids (only 'yes' means yes).
filtered_bam_root=$5 # root name for output bam (e.g. "out" will create "out.bam" and "out_flagstat.txt") 

unfiltered_bam_root=${unfiltered_bam%.bam}
echo "-- Filtered alignments file will be: '${filtered_bam_root}.bam'"

echo "-- Sort bam by name."
set -x
samtools sort -@ $ncpus -m 4G -n -O sam -T sorted $unfiltered_bam > sorted.sam
set +x
echo "-- Handle UMI flagging and errors with 'filter_reads.py'."
# NOTE script written for python3 works just as well for python2.7 as long as pysam works
marked_root="${unfiltered_bam_root}_marked"
set -x
python2.7 ./filter_reads.py --min_mapq $map_thresh sorted.sam flagged_presorted.sam
set +x
#samtools view -bS flagged_presorted.sam > flagged.bam
echo "-- Sort bam by location."
set -x
samtools sort -@ $ncpus -m 4G -O bam -T flagged flagged_presorted.sam > flagged.bam
rm *.sam
set +x

echo "-- Add mate cigar information..."
set -x
time java -jar ./picard.jar RevertOriginalBaseQualitiesAndAddMateCigar \
    INPUT=flagged.bam OUTPUT=cigar.bam VALIDATION_STRINGENCY=SILENT \
    RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0
set +x

if [ "$umi" == "yes" ] || [ "$umi" == "y" ] || [ "$umi" == "true" ] || [ "$umi" == "t" ] || [ "$umi" == "umi" ]; then

    # MarkDuplicates on UMI is non-comparable (to non-UMI) and confusing.
	
    echo "-- Running picard mark duplicates on UMI..."
    # From:
	# stampipes/makefiles/picard/dups_cigarumi.mk
    set -x
    time java -jar ./picard.jar UmiAwareMarkDuplicatesWithMateCigar \
        INPUT=cigar.bam OUTPUT=marked.bam METRICS_FILE=${filtered_bam_root}_dup_qc.txt \
        UMI_TAG_NAME=XD ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'    
	set +x
		
	echo "-- ------------- umi_dups"
	cat ${filtered_bam_root}_umi_dups.txt
	echo "-- -------------"
    ls -l
    # Richard Sandstrom: UMI flags: 512 + 1024 = 1536 (both 8 and 4 are incorporated into 512 by the script. 1024 is set unambiguously by the UMI eval portion of the script)
    echo "-- UMI filtering will be performed."
    filter_flags=`expr 512 + 1024`

else
    echo "-- Running picard mark duplicates on non-UMI..."
    # From:
	# stampipes/makefiles/picard/dups_cigar.mk
    set -x
    time java -jar ./picard.jar MarkDuplicatesWithMateCigar \
        INPUT=cigar.bam OUTPUT=marked.bam METRICS_FILE=${filtered_bam_root}_dup_qc.txt \
	    ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
        READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
    set +x

    # Richard Sandstrom: non-UMI flags: 512 only (again, 8 and 4 are both criteria to set 512.  we don't filter dups for non UMI reads by convention).
    filter_flags=512
fi

if [ -e ${filtered_bam_root}_dup_qc.txt ]; then
	echo "-- ------------- picard MarkDuplicates"
	cat ${filtered_bam_root}_dup_qc.txt
	echo "-- -------------"
fi

echo "-- Filter bam and threshold..."
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
set -x
samtools view -F $filter_flags -q ${map_thresh} -b marked.bam > ${filtered_bam_root}.bam
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

