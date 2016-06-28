#!/bin/bash -e

if [ $# -ne 8 ]; then
    echo "usage v1: dnase_align_bwa_pe.sh <index.tgz> <reads1.fq.gz> <reads2.fq.gz> <barcode> <umi> <adapters.tsv> <ncpus> <bam_root>"
    echo "Align paired-end reads with bwa.  Is independent of DX and encodeD."
    echo "Requires bwa, edwBamStats, and samtools on path."
    exit -1; 
fi
index_tgz=$1  # BWA index archive including <ref_id>.fa (e.g. GRCh38.fa) and bwa index.
reads1_fq_gz=$2  # fastq of of paired-end read1, which will be trimmed resulting in "read1_trimmed.fq.gz"
reads2_fq_gz=$3  # fastq of of paired-end read2, which will be trimmed resulting in "read2_trimmed.fq.gz"
barcode=$4       # Barcode used in generating fastqs
umi=$5           # Whether reads in bam contain UMI ids (only 'yes' means yes).
adapters_tsv=$6  # Adapter file that contains <barcode>_P5 and <barcode>_P7 sequences that should be trimmed
ncpus=$7         # Number of cpus available.
bam_root=$8      # root name for output bam (e.g. "out" will create "out.bam" and "out_flagstat.txt") 

echo "-- Expect to create '${bam_root}.bam'"

if [ "$umi" == "yes" ] || [ "$umi" == "y" ] || [ "$umi" == "true" ] || [ "$umi" == "t" ] || [ "$umi" == "umi" ]; then
    echo "-- UMI marking fastqs..."
    set -x
    python2.7 /usr/bin/fastq_umi_add.py $reads1_fq_gz R1_umi.fq.gz
    python2.7 /usr/bin/fastq_umi_add.py $reads2_fq_gz R2_umi.fq.gz
    set +x
    reads1_fq_gz=R1_umi.fq.gz
    reads2_fq_gz=R2_umi.fq.gz
fi

echo "-- Adapter trimming..."
set -x
trim-adapters-illumina -f $adapters_tsv -1 ${barcode}_P5 -2 ${barcode}_P7 $reads1_fq_gz $reads2_fq_gz \
																R1_trimmed.fq.gz R2_trimmed.fq.gz 2> ${bam_root}_trim_stats.txt
set +x
reads1_fq_gz=R1_trimmed.fq.gz
reads2_fq_gz=R2_trimmed.fq.gz

# --- Debug only
echo "-- Trimming stats..."
set -x
cat ${bam_root}_trim_stats.txt 
set +x
echo "--------------------"

echo "-- Uncompress index archive..."
set -x
tar zxvf $index_tgz
set +x
ref_fa=`ls *.fa`
ref_id=${ref_fa%.fa}

echo "-- Aligning with bwa to ${ref_id}..."
set -x
bwa aln -q 5 -l 32 -k 2 -t $ncpus $ref_id $reads1_fq_gz > tmp_1.sai
bwa aln -q 5 -l 32 -k 2 -t $ncpus $ref_id $reads2_fq_gz > tmp_2.sai
bwa sampe $ref_id tmp_1.sai tmp_2.sai $reads1_fq_gz $reads2_fq_gz \
                | samtools view -Shu - | samtools sort -@ $ncpus -m 5G -f - tmp.sam
samtools view -hb tmp.sam > ${bam_root}.bam
samtools index ${bam_root}.bam
#samtools view -H ${bam_root}.bam
set +x

#echo "-- Clean-up..."
#set -x
#rm -f tmp.sai tmp.bam
#rm -f ${ref_id}.*
#set +x

echo "-- Collect bam stats..."
set -x
samtools flagstat ${bam_root}.bam > ${bam_root}_flagstat.txt
edwBamStats ${bam_root}.bam ${bam_root}_edwBamStats.txt
set +x

echo "-- The results..."
ls -l ${bam_root}*
df -k .

