#!/bin/bash -e

if [ $# -ne 8 ]; then
    echo "usage v1: dnase_qc_bam.sh <filtered.bam> <sample_size> <ncpus> <pe_or_se> <assembly> <mappable.tgz> <read_length> <hotspot_dir>"
    echo "Evaluates filtered aligned reads for DNase using HotSpot1 and other tools.  Is independent of DX and encodeD."
    echo "Requires edwBamFilter, edwBamStats, samtools, Rscript, phantompeakqualtools, caTools, snow, spp, gawk,"
    echo "         hotspot, hotspot.py (GCAP), bedops and bedtools on path."
    exit -1; 
fi
filtered_bam=$1 # filtered bam file.
sample_size=$2  # number of sample reads to evaluate (e.g. 15000000)
ncpus=$3        # Number of cpus available
pe_or_se=$4     # Either "pe" for paired-end or "se" for single-end
assembly=$5     # Genome assembly matching chrom.sizes and bam.
mappable_tgz=$6 # Archive with mappable_target.starch, center_sites.starch and chrom_sizes.bed created by dnase_index.bwa.sh.
read_length=$7  # Length of reads to match hotspot Kmer size (must be 32, 36, 40, 50, 58, 72, 76, or 100).
hotspot_dir=$8  # Directory where hotpot is installed 
if [ "$pe_or_se" != "pe" ] && [ "$pe_or_se" != "se" ]; then
    echo "-- ERROR: Declare either 'pe' for paied-end or 'se' for single-end."
    exit 1
fi
if [ "$read_length" != "32" ] && [ "$read_length" != "36" ] && [ "$read_length" != "40" ] && [ "$read_length" != "50" ] \
&& [ "$read_length" != "58" ] && [ "$read_length" != "72" ] && [ "$read_length" != "76" ] && [ "$read_length" != "100" ]; then
    echo "* ERROR: Read length ($read_length) must be one of 32, 36, 40, 50, 58, 72, 76, or 100."
    exit 1
fi
if [ "$assembly" == "GRCh38" ]; then
    # HotSpot1 uses UCSC assembly names
    assembly="hg38"
fi 

bam_input_root=${filtered_bam%.bam}
bam_no_chrM_root="${bam_input_root}_no_chrM"
bam_sample_root="${bam_input_root}_${sample_size}_sample"
echo "-- Sampled alignments file will be: '${bam_sample_root}.bam'"

echo "-- Filter out chrM..."
set -x
edwBamFilter -sponge -chrom=chrM $filtered_bam ${bam_no_chrM_root}.bam  ## qc based on bam without chrm
set +x
if [ "$pe_or_se" == "pe" ]; then
    echo "-- Sorting by name for paired-end data..."
    # Note the sort by name which is needed for proper pe sampling
    set -x
    samtools sort -@ $ncpus -m 6G -n -f ${bam_no_chrM_root}.bam ${bam_no_chrM_root}_byname.sam ## for pbc usage
    samtools view -hb ${bam_no_chrM_root}_byname.sam > ${bam_no_chrM_root}_byname.bam
    rm *.sam
    set +x
    bam_no_chrM_root="${bam_no_chrM_root}_byname"
fi
set -x
samtools index ${bam_no_chrM_root}.bam
set +x

echo "-- Generating stats on $sample_size reads..."
set -x
edwBamStats -sampleBamSize=$sample_size -u4mSize=$sample_size -sampleBam=${bam_sample_root}.bam \
                      ${bam_no_chrM_root}.bam ${bam_no_chrM_root}_sampling_edwBamStats.txt
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
bedtools bamtobed -i ${bam_sample_root}.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c \
    | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
    > ${bam_sample_root}_pbc.txt
set +x

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

echo "-- Moving hotspot1 mappable file to expected name and place..."
mappable=${assembly}.K${read_length}.mappable_only
satellites=Satellite.${assembly}
if [ "$read_length" == "36" ]; then
    set -x
    unstarch mappable_target.starch > ${hotspot_dir}/hotspot-distr/data/${mappable}.bed
    # Note that the blacklist is already removed from mappable_targets so will be a noop but hotspot1 doesn't know that.
    if [ -e ${assembly}.blacklist.bed ]; then
        set -x
        cp ${assembly}.blacklist.bed ${hotspot_dir}/hotspot-distr/data/${satellites}.bed
        set +x
    else    # hg19 has wgEncodeDacMapabilityConsensusExcludable.bed instead
        set -x
        touch ${hotspot_dir}/hotspot-distr/data/${satellites}.bed
        set +x
    fi
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
hotspot.py -o ${hotspot_dir}/hotspot-distr/ ${bam_sample_root}.bam $assembly DNase-seq $read_length tmp out

cp tmp/${bam_sample_root}.spot.out ${bam_sample_root}_hotspot1_qc.txt
set +x

echo "-- The results..."
ls -l ${bam_sample_root}*
df -k .

