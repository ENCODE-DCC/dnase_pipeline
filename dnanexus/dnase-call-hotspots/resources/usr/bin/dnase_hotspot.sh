#!/bin/bash -e

if [ $# -lt 7 ] || [ $# -gt 8 ]; then
    echo "usage v1: dnase_hotspot.sh <aligned.bam> <chrom.sizes> <mappable.tgz> <hotspot_root> <peaks_root> <density_root> <include_minor_chroms> [<allcalls_root>]"
    echo "Calls peaks hotspot2.  Is independent of DX and encodeD."
    echo "Requires hotspot2 (hotspot2,hotspot2.sh,cutcounts.bash,density-peaks.bash,hsmerge.sh),"
    echo "         bedops (bam2bed,bedmap,sort-bed,starch,unstarch), modwt, mawk, samtools, bedToBigBed,"
    echo "         bedGraphToBigWig, and pigz on path."
    exit -1;
fi
bam_file=$1      # Bam file on which hotspot will be run.
chrom_sizes=$2   # Chrom sizes file that matches the reference genome which bam reads were aligned to.
mappable_tgz=$3  # Archive with mappable_target.starch, center_sites.starch and chrom_sizes.bed created by dnase_index.bwa.sh.
hotspot_root=$4  # Put hotspot results into ${hotspot_root}.bed.gz, ${hotspot_root}.bb, ${hotspot_root}_count.txt, and ${hotspot_root}_SPOT.txt
peaks_root=$5    # Put peak results into ${peaks_root}.bed.gz, ${peaks_root}.bb, and ${peaks_root}_count.txt
density_root=$6  # Put density results into ${density_root}.bw
minor_chroms=$7  # if 'true' call hotspots on chrM and scaffolds as well as the default 1-22+XY chromosomes.
echo "-- output: '${hotspot_root}.bed.gz/.bb/_count.txt/_SPOT.txt', '${peaks_root}.bed.gz/.bb/_count.txt', '${density_root}.bw'"
if [ $# -eq 8 ]; then
    allcalls_root=$8 # (optional) Put all_calls results into ${allcalls_root}.bed.gz and ${allcalls_root}_count.txt
    echo "           and '${allcalls_root}.bed.gz/_count.txt'"
fi

echo "-- Extracting mappable regions..."
set -x
tar -xzf $mappable_tgz
set +x

if [ "$minor_chroms" != "true" ]; then
    echo "-- Filtering bam and mappables to only 22+XY..."
    chroms='chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY'
    chroms_jailed=`echo $chroms | sed 's/ /|/g'`
    set -x
    mv $bam_file input.bam
    mv center_sites.starch centers_input.starch
    mv mappable_target.starch mappable_input.starch
    mv chrom_sizes.bed chrom_sizes_input.bed
    samtools index input.bam
    samtools view -b input.bam $chroms > $bam_file
    unstarch centers_input.starch | grep -w -E $chroms_jailed | starch - > center_sites.starch
    unstarch mappable_input.starch | grep -w -E $chroms_jailed | starch - > mappable_target.starch
    grep -w -E $chroms_jailed chrom_sizes_input.bed > chrom_sizes.bed
    rm input.bam* centers_input.starch mappable_input.starch chrom_sizes_input.bed
    set +x
fi

echo "-- Running hotspot2.sh..."
set -x
hotspot2.sh -c chrom_sizes.bed -C center_sites.starch -M mappable_target.starch $bam_file out/
set +x

echo "-- hotspot2.sh out/..."
ls -l out

echo "-- Converting hotspots to bed.gz and bigBed..."
set -x
mv out/*.hotspots.fdr*.starch ${hotspot_root}.starch
unstarch --elements ${hotspot_root}.starch > ${hotspot_root}_count.txt
unstarch ${hotspot_root}.starch > ${hotspot_root}.bed
bedToBigBed -as=/usr/bin/broadPeak.as -type=bed6+3 ${hotspot_root}.bed $chrom_sizes ${hotspot_root}.bb
#grep "^chr" ${hotspot_root}.bed | wc -l > ${hotspot_root}_count.txt
pigz ${hotspot_root}.bed
set +x

echo "-- Converting narrowpeaks to bed.gz and bigBed..."
set -x
mv out/*.peaks.narrowpeaks.starch ${peaks_root}.starch
unstarch --elements ${peaks_root}.starch > ${peaks_root}_count.txt
unstarch ${peaks_root}.starch > ${peaks_root}.bed
bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${peaks_root}.bed $chrom_sizes ${peaks_root}.bb
#grep "^chr" ${peaks_root}.bed | wc -l > ${peaks_root}_count.txt
pigz ${peaks_root}.bed
set +x

echo "-- Saving density bigWig..."
set -x
mv out/*.density.bw ${density_root}.bw
mv out/*.density.starch ${density_root}.starch
set +x

echo "-- Saving SPOT..."
set -x
mv out/*.SPOT.txt ${hotspot_root}_SPOT.txt
set +x
echo "-- ----- SPOT score: -----"
cat ${hotspot_root}_SPOT.txt

if [ $# -eq 8 ]; then
    echo "-- Counting allcalls..."
    set -x
    mv out/*.allcalls.starch ${allcalls_root}.starch
    unstarch --elements ${allcalls_root}.starch > ${allcalls_root}_count.txt
    unstarch ${allcalls_root}.starch > ${allcalls_root}.bed
    #grep "^chr" ${allcalls_root}.bed | wc -l > ${allcalls_root}_count.txt
    pigz ${allcalls_root}.bed
    set +x
fi

echo "-- The results..."
ls -l ${hotspot_root}*
ls -l ${peaks_root}*
ls -l ${density_root}*
if [ $# -eq 8 ]; then
    ls -l ${allcalls_root}*
fi
df -k .


