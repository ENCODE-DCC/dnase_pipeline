#!/bin/bash -e

if [ $# -ne 6 ]; then
    echo "usage v1: dnase_qc_hotspot.sh <full.bam> <chrom.sizes> <genome> <sample_size> <read_length> <hotspot_dir>"
    echo "Performs hotspot on a sample of a bam for QC results.  Is independent of DX and encodeD."
    echo "Requires hotspot, hotspot.py (GCAP), samtools, edwBamStats, "
    echo "         bedops (bedmap,sort-bed,starch,starchcat,unstarch), and bedtools (bamToBed,shuffleBed) on path."
    exit -1; 
fi
full_bam=$1     # Bam file from which sample will be taken.
chrom_sizes=$2  # Chrom sizes file that matches the reference genome which bam reads were aligned to.
genome=$3       # Genome assembly matching chrom.sizes and bam.
sample_size=$4  # number of reads to sample from bam for input into hotspot.
read_length=$5  # Length of reads to match hotspot Kmer size (must be 32, 36, 40, 50, 58, 72, 76, or 100).
hotspot_dir=$6  # Directory where hotpot is installed 
#                 (e.g. if "/usr/local/hotspot" then "/usr/local/hotspot/hotspot-distr/hotspot-deploy/bin" should be in path)

if [ "$read_length" != "32" ] && [ "$read_length" != "36" ] && [ "$read_length" != "40" ] && [ "$read_length" != "50" ] \
&& [ "$read_length" != "58" ] && [ "$read_length" != "72" ] && [ "$read_length" != "76" ] && [ "$read_length" != "100" ]; then
    echo "* ERROR: Read length ($read_length) must be one of 32, 36, 40, 50, 58, 72, 76, or 100."
    exit 1
fi

full_root=${full_bam%.bam}
sample_root="${full_root}_sample_${sample_size}"
echo "-- sampled bam will be named: '${sample_root}.bam'"

echo "-- Creating chromInfo.bed from ${chrom_sizes}..."
grep -w chr[1-9] $chrom_sizes > min_chrom.sizes
grep -w chr[1-2][0-9] $chrom_sizes >> min_chrom.sizes
grep -w chr[X,Y] $chrom_sizes >> min_chrom.sizes
# sort-bed is important!
cat min_chrom.sizes | awk '{printf "%s\t0\t%s\t%s\n",$1,$2,$1}' | sort-bed - > ${genome}.chromInfo.bed

if [ ! -f ${full_bam}.bai ]; then
    echo "-- Indexing bam..."
    set -x
    samtools index $full_bam
    set +x
fi
    
echo "-- Sampling bam..."
set -x
edwBamStats -sampleBamSize=$sample_size -u4mSize=$sample_size -sampleBam=${sample_root}.bam $full_bam \
                                                                                  ${full_root}_sampling_edwBamStats.txt
edwBamStats ${sample_root}.bam ${sample_root}_edwBamStats.txt
samtools index ${sample_root}.bam
set +x
    
echo "-- Running hotspot.py..."
set -x
mkdir tmp
mkdir out
cp ${genome}.chromInfo.bed ${hotspot_dir}/hotspot-distr/data/
# TODO:
# May also need to do something about "Satellite.${genome}.bed"
#cp /usr/bin/Satellite.${genome}.bed ${hotspot_dir}/hotspot-distr/data   # hg19 version already there!
mappable=${genome}.K${read_length}.mappable_only
wget http://www.uwencode.org/proj/hotspot/${mappable}.bed -O ${hotspot_dir}/hotspot-distr/data/${mappable}.bed
hotspot.py -o ${hotspot_dir}/hotspot-distr/ ${sample_root}.bam $genome DNase-seq $read_length tmp out
# hanging on: /home/dnanexus/hotspot/hotspot-distr/pipeline-scripts/run_generate_random_lib

cp tmp/${sample_root}.spot.out ${sample_root}_hotspot_qc.txt
set +x

echo "-- The results..."
ls -l ${sample_root}*
df -k .

