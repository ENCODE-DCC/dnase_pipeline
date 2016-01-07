#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: dnase_align_bwa_pe.sh <index.tgz> <reads1.fq.gz> <reads2.fq.gz> <ncpus> <bam_root>"
    echo "Align paired-end reads with bwa.  Is independent of DX and encodeD."
    echo "Requires bwa, edwBamStats, and samtools on path."
    exit -1; 
fi
index_tgz=$1  # BWA index archive including <ref_id>.fa (e.g. GRCh38.fa) and bwa index.
reads1_fq_gz=$2  # fastq of of paired-end read1, which will be trimmed resulting in "read1_trimmed.fq.gz"
reads2_fq_gz=$3  # fastq of of paired-end read2, which will be trimmed resulting in "read2_trimmed.fq.gz"
ncpus=$4      # Number of cpus available.
bam_root="${5}_pe_bwa"   # root name for output bam (e.g. "out" will create "out_pe_bwa.bam" and "out_pe_bwa_flagstat.txt") 

echo "-- Expect to create '${bam_root}.bam'"

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

