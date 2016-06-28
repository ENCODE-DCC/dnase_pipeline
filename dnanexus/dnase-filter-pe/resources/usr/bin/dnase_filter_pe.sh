#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: dnase_filter_pe.sh <unfiltered.bam> <map_threshold> <ncpus> <umi> <filtered_bam_root>"
    echo "Filters paired-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires: filter_reads.py, mark_duplicates.mk and umi_sort_sam_annotate.awk in /usr/bin; and samtools on path."
    echo "          Also needs pysam, java and picard 1.92."
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
samtools sort -@ $ncpus -m 6G -n -O sam -T sorted $unfiltered_bam > sorted.sam
set +x
echo "-- Handle UMI flagging and errors with 'filter_reads.py'."
# NOTE script written for python3 works just as well for python2.7 as long as pysam works
marked_root="${unfiltered_bam_root}_marked"
set -x
python2.7 /usr/bin/filter_reads.py --min_mapq $map_thresh sorted.sam flagged_presorted.sam
#samtools view -bS flagged_presorted.sam > flagged.bam
samtools sort -@ $ncpus -m 6G -O bam -T flagged flagged_presorted.sam > flagged.bam
rm *.sam
set +x
    # sam to bam?

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

    # MarkDuplicates on UMI is non-comparable (to non-UMI) and confusing.
	
    echo "-- UMI filtering will be performed."
    set -x
    mkdir tmp
    make -f /usr/bin/mark_umi_dups.mk THREADS=$ncpus MAX_MEM=6 TMPDIR=tmp INPUT_BAM_FILE=flagged.bam OUTPUT_BAM_FILE=marked.bam
    set +x

    echo "-- Counting UMI marked duplicates..."
    set -x
	all=`samtools view -c -F 512 marked.bam`
	dups=`samtools view -c -F 512 -f 1024 marked.bam`
	pct_dups=`awk '{printf("%.6f\n",$1/$2)}' <<<" ${dups} ${all}"`
    set +x
	# Match picard_MarkDuplicates formatting so that parsing works the same.
	echo -e "From\tReads_Examined\tRead_Duplicates\tUMI_Read_Duplicates\tPercent_Duplication" > ${filtered_bam_root}_umi_dups.txt
	echo -e "samtools view -c\t${all}\t${dups}\t${dups}\t${pct_dups}"                        >> ${filtered_bam_root}_umi_dups.txt

	echo "-- ------------- umi_dups"
	cat ${filtered_bam_root}_umi_dups.txt
	echo "-- -------------"
    ls -l
    # Richard Sandstrom: UMI flags: 512 + 1024 = 1536 (both 8 and 4 are incorporated into 512 by the script. 1024 is set unambiguously by the UMI eval portion of the script)
    filter_flags=`expr 512 + 1024`

else
    echo "-- Running picard mark duplicates on non-UMI..."
    set -x
    time java -Xmx4G -jar /picard/MarkDuplicates.jar INPUT=flagged.bam OUTPUT=marked.bam \
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
fi

echo "-- Filter bam and threshold..."
set -x
samtools view -F $filter_flags -q ${map_thresh} -b marked.bam > ${filtered_bam_root}.bam
set +x
#set -x
#samtools view -F $filter_flags -f 2 -q ${map_thresh} -u marked.bam | \
#        samtools sort -@ $ncpus -m 6G -n -f - ${filtered_bam_root}_tmp.sam
#samtools view -hb ${filtered_bam_root}_tmp.sam > ${filtered_bam_root}_tmp.bam
#samtools fixmate -r ${filtered_bam_root}_tmp.bam -O sam - | \
#        samtools view -F $filter_flags -f 2 -u - | \
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

