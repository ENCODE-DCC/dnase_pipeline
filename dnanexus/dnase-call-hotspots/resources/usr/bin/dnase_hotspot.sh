#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: dnase_hotspot.sh <aligned.bam> <chrom.sizes> <genome> <read_length> <hotspot_dir>"
    echo "Calls peaks and regions with hotspot.  Is independent of DX and encodeD."
    echo "Requires hotspot, hotspot.py (GCAP), samtools, bedToBigBed, bedGraphToBigWig, bedGraphPack, edwBamStats"
    echo "         bedops (bedmap,sort-bed,starch,starchcat,unstarch), and bedtools (bamToBed,intersectBed,shuffleBed) on path."
    exit -1; 
fi
bam_file=$1     # Bam file on which hotspot will be run.
chrom_sizes=$2  # Chrom sizes file that matches the reference genome which bam reads were aligned to.
genome=$3       # Genome assembly matching chrom.sizes and bam.
read_length=$4  # Length of reads to match hotspot Kmer size (must be 32, 36, 40, 50, 58, 72, 76, or 100).
hotspot_dir=$5  # Directory where hotpot is installed 
#                 (e.g. if "/usr/local/hotspot" then "/usr/local/hotspot/hotspot-distr/hotspot-deploy/bin" should be in path)

if [ "$read_length" != "32" ] && [ "$read_length" != "36" ] && [ "$read_length" != "40" ] && [ "$read_length" != "50" ] \
&& [ "$read_length" != "58" ] && [ "$read_length" != "72" ] && [ "$read_length" != "76" ] && [ "$read_length" != "100" ]; then
    echo "* ERROR: Read length ($read_length) must be one of 32, 36, 40, 50, 58, 72, 76, or 100."
    exit 1
fi

bam_root=${bam_file%.bam}
narrowPeak_root="${bam_root}_narrowPeak_hotspot"
broadPeak_root="${bam_root}_broadPeak_hotspot"
signal_root="${bam_root}_signal_hotspot"
echo "-- output: '${narrowPeak_root}.bed/.bb', '${broadPeak_root}.bed/.bb', '${signal_root}.bw'"

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
    
echo "-- Running hotspot.py..."
set -x
mkdir tmp
mkdir out
cp ${genome}.chromInfo.bed ${hotspot_dir}/hotspot-distr/data/
# May also need to do something about "Satellite.${genome}.bed"
#cp /usr/bin/Satellite.${genome}.bed ${hotspot_dir}/hotspot-distr/data   # hg19 version already there!
mappable=${genome}.K${read_length}.mappable_only
wget http://www.uwencode.org/proj/hotspot/${mappable}.bed -O ${hotspot_dir}/hotspot-distr/data/${mappable}.bed
hotspot.py ${hotspot_dir}/hotspot-distr/ ${bam_root}.bam $genome DNase-seq $read_length tmp out
cp tmp/${bam_root}.spot.out ${bam_root}_hotspot_out.txt
set +x

echo "-- Generating narrowPeaks..."
set -x
# Convert output into ENCODE formats,
# Convert narrowPeaks.bed and several stray columns to narrowPeaks.bigBed in $4
paste out/narrowPeaks.bed out/narrowPeaks.dens out/narrowPeaks.pval | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "p" NR, 0, ".", $4, $5, -1, -1}' - | \
    sort -k1,1 -k2,2n - > ${narrowPeak_root}.bed
bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${narrowPeak_root}.bed $chrom_sizes ${narrowPeak_root}.bb
wc -l ${narrowPeak_root}.bed > ${narrowPeak_root}_qc.txt
set +x

echo "-- Generating broadPeaks..."
set -x
# Convert broadPeaks.bed and broadPeaks.pVal to broadPeaks.bigBed in $5
paste out/broadPeaks.bed out/broadPeaks.pval | \
    awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "hot" NR, 0, ".", $5, $6, -1}' - | \
    sort -k1,1 -k2,2n - > ${broadPeak_root}.bed
bedToBigBed -as=/usr/bin/broadPeak.as -type=bed6+3 ${broadPeak_root}.bed $chrom_sizes ${broadPeak_root}.bb
wc -l ${broadPeak_root}.bed > ${broadPeak_root}_qc.txt
set +x

echo "-- Generating signal..."
set -x
# Convert starched bedGraph to mappable-only bigWig in $6
unstarch out/density.bed.starch > tmp/tmp.bed
intersectBed -a tmp/tmp.bed -b ${hotspot_dir}/hotspot-distr/data/${genome}.K${read_length}.mappable_only.bed -f 1.00 | \
    cut -f 1,2,3,5 | bedGraphPack stdin ${signal_root}.bedGraph
bedGraphToBigWig ${signal_root}.bedGraph $chrom_sizes ${signal_root}.bw
set +x

echo "-- The results..."
ls -l ${bam_root}_hotspot_out.txt ${narrowPeak_root}* ${broadPeak_root}* ${signal_root}*
df -k .

