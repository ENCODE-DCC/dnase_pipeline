#!/bin/bash -e

if [ $# -ne 5 ]; then
    echo "usage v1: dnase_qc_hotspot1.sh <file.bam> <assembly> <mappable.tgz> <read_length> <hotspot_dir>"
    echo "Performs hotspot1 on a (sampled) bam for QC results.  Is independent of DX and encodeD."
    echo "Requires hotspot, hotspot.py (GCAP), samtools, "
    echo "         bedops (bedmap,sort-bed,starch,starchcat,unstarch), and bedtools (bamToBed,shuffleBed) on path."
    exit -1; 
fi
bam_file=$1     # Bam file from which sample will be taken.
assembly=$2     # Genome assembly matching chrom.sizes and bam.
mappable_tgz=$3 # Archive with mappable_target.starch, center_sites.starch and chrom_sizes.bed created by dnase_index.bwa.sh.
read_length=$4  # Length of reads to match hotspot Kmer size (must be 32, 36, 40, 50, 58, 72, 76, or 100).
hotspot_dir=$5  # Directory where hotpot is installed 
#                 (e.g. if "/usr/local/hotspot" then "/usr/local/hotspot/hotspot-distr/hotspot-deploy/bin" should be in path)

if [ "$read_length" != "32" ] && [ "$read_length" != "36" ] && [ "$read_length" != "40" ] && [ "$read_length" != "50" ] \
&& [ "$read_length" != "58" ] && [ "$read_length" != "72" ] && [ "$read_length" != "76" ] && [ "$read_length" != "100" ]; then
    echo "* ERROR: Read length ($read_length) must be one of 32, 36, 40, 50, 58, 72, 76, or 100."
    exit 1
fi
if [ "$assembly" == "GRCh38" ]; then
    # HotSpot1 uses UCSC assembly names
    assembly="hg38"
fi 

bam_root=${bam_file%.bam}
echo "-- qc results will be named: '${bam_root}_hotspot1_qc.txt'"

echo "-- Opening archive of mappable regions..."
set -x
tar -xzf $mappable_tgz
set +x
# Expect: chrom_sizes.bed, center_sites.starch, mappable_target.starch, GRCh38_no_alts.K36.mappable_only.starch, GRCh38.blacklist.bed

echo "-- Creating chromInfo.bed from chrom_sizes.bed..."
# Reducing to min chroms appears to be necessary for hotspot1
grep -w chr[1-9] chrom_sizes.bed > min_chrom.sizes
grep -w chr[1-2][0-9] chrom_sizes.bed >> min_chrom.sizes
grep -w chr[X,Y] chrom_sizes.bed >> min_chrom.sizes
## sort-bed is important!
cat min_chrom.sizes | awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,$1}' | sort-bed - > ${assembly}.chromInfo.bed

if [ ! -f ${bam_file}.bai ]; then
    echo "-- Indexing bam..."
    set -x
    samtools index $bam_file
    set +x
fi
    
echo "-- Moving hotspot1 mappable file to expected name and place..."
mappable=${assembly}.K${read_length}.mappable_only
satellites=Satellite.${assembly}
if [ "$read_length" == "36" ]; then
    set -x
    unstarch mappable_target.starch > ${hotspot_dir}/hotspot-distr/data/${mappable}.bed
    # Note that the blacklist is already removed from mappable_targets so will be a noop but hotspot1 doesn't know that.
    cp *.blacklist.bed ${hotspot_dir}/hotspot-distr/data/${satellites}.bed
    set +x
else
    # Not a good idea: probably only hg19 will work
    set -x
    wget http://www.uwencode.org/proj/hotspot/${mappable}.bed -O ${hotspot_dir}/hotspot-distr/data/${mappable}.bed
    wget http://www.uwencode.org/proj/hotspot/${satellites}.bed -O ${hotspot_dir}/hotspot-distr/data/${satellites}.bed
    set +x
fi

echo "-- Running hotspot.py..."
set -x
mkdir tmp
mkdir out
cp ${assembly}.chromInfo.bed ${hotspot_dir}/hotspot-distr/data/
hotspot.py -o ${hotspot_dir}/hotspot-distr/ ${bam_file} $assembly DNase-seq $read_length tmp out
# hanging on: /home/dnanexus/hotspot/hotspot-distr/pipeline-scripts/run_generate_random_lib

cp tmp/${bam_root}.spot.out ${bam_root}_hotspot1_qc.txt
set +x

echo "-- The results..."
ls -l ${bam_root}*
df -k .

