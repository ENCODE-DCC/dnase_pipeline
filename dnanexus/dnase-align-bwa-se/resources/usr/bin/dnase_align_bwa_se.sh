#!/bin/bash -e

if [ $# -ne 4 ]; then
    echo "usage v1: dnase_align_bwa_se.sh <index.tgz> <reads.fq.gz> <ncpus> <bam_root>"
    echo "Align single-end reads with bwa.  Is independent of DX and encodeD."
    echo "Requires bwa, edwBamStats, and samtools on path."
    exit -1; 
fi
index_tgz=$1   # BWA index archive including <ref_id>.fa (e.g. GRCh38.fa) and bwa index.
reads_fq_gz=$2 # fastq of of single-end reads, which will be trimmed resulting in "read1_trimmed.fq.gz"
ncpus=$3       # Number of cpus.
bam_root=$4    # root name for output bam (e.g. "out" will create "out.bam" and "out_flagstat.txt") 

echo "-- Expect to create '${bam_root}.bam'"

echo "-- Uncompress index archive..."
set -x
tar zxvf $index_tgz
set +x
ref_fa=`ls *.fa`
ref_id=${ref_fa%.fa}

echo "-- Aligning with bwa..."
set -x
bwa aln -Y -l 32 -n 0.04 -k 2 -t $ncpus $ref_id $reads_fq_gz > tmp.sai
bwa samse -n 10 ${ref_id} tmp.sai $reads_fq_gz | samtools view -Shu - | samtools sort -@ $ncpus -m 5G -f - tmp.sam
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

