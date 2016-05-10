#!/bin/bash -e

if [ $# -ne 3 ]; then
    echo "usage v1: dnase_eval_pe.sh <filtered.bam> <sample_size> <ncpus>"
    echo "Evaluates filtered paired-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires edwBamFilter,edwBamStats,samtools,Rscript,phantompeakqualtools,caTools,snow,spp,gawk,bedtools on path."
    exit -1; 
fi
filtered_bam=$1  # filtered bam file.
sample_size=$2   # number of sample reads to evaluate (e.g. 15000000)
ncpus=$3         # Number of cpus available

bam_input_root=${filtered_bam%.bam}
bam_no_chrM_root="${bam_input_root}_no_chrM"
bam_sample_root="${bam_input_root}_${sample_size}_sample"
echo "-- Sampled alignments file will be: '${bam_sample_root}.bam'"

echo "-- Filter out chrM..."
# Note the sort by name which is needed for proper pe sampling
set -x
edwBamFilter -sponge -chrom=chrM $filtered_bam ${bam_no_chrM_root}.bam  ## qc based on bam without chrm
samtools sort -@ $ncpus -m 6G -n -f ${bam_no_chrM_root}.bam ${bam_no_chrM_root}_byname.sam ## for pbc usage
samtools view -hb ${bam_no_chrM_root}_byname.sam > ${bam_no_chrM_root}_byname.bam
samtools index ${bam_no_chrM_root}_byname.bam
rm *.sam
set +x

echo "-- Generating stats on $sample_size reads..."
set -x
edwBamStats -sampleBamSize=$sample_size -u4mSize=$sample_size -sampleBam=${bam_sample_root}.bam \
                      ${bam_no_chrM_root}_byname.bam ${bam_no_chrM_root}_sampling_edwBamStats.txt
edwBamStats ${bam_sample_root}.bam ${bam_sample_root}_edwBamStats.txt
samtools index ${bam_sample_root}.bam
set +x

echo "-- Running spp..."
set -x
# awk didn't work, so use gawk and pre-create the tagAlign 
samtools view -F 0x0204 -o - ${bam_sample_root}.bam | \
   gawk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' \
   | gzip -c > ${bam_sample_root}.tagAlign.gz
Rscript phantompeakqualtools/run_spp.R -x=-500:-1 -s=-500:5:1500 -rf -c=${bam_sample_root}.tagAlign.gz \
                                    -out=${bam_sample_root}_spp.txt > ${bam_sample_root}_spp_out.txt
touch ${bam_sample_root}_spp.txt
set +x

echo "-- Running pbc..."
# Seth interprets:
# TotalReadPairs \t DistinctReadPairs \t OneReadPair \t TwoReadPairs \t NRF=Distinct/Total \t PBC1=OnePair/Distinct \t PBC2=OnePair/TwoPair
# Tim reinterprets:
# Sampled Reads \t Distinct Locations Mapped \t Single-read Locations \t Multi-read Locations \t NRF (Non-Redundant Fraction)=Distinct Locations/Sample Reads \t PBC1=Single-read Locations/Distinct Locations \t PBC2=Single-read Locations/Multi-read Locations
set -x
bedtools bamtobed -i ${bam_sample_root}.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
    | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > ${bam_sample_root}_pbc.txt
set +x

echo "-- The results..."
ls -l ${bam_sample_root}*
df -k .

