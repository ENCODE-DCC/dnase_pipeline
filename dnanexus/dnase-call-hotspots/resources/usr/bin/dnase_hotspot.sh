#!/bin/bash -e

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "usage v1: dnase_hotspot.sh <aligned.bam> <chrom.sizes> <blacklist.bed>"
    echo "Calls peaks hotspot2.  Is independent of DX and encodeD."
    echo "Requires hotspot2 (hotspot2,hotspot2.sh,tallyCountsInSmallWindows,cutcounts.bash,density-peaks.bash), "
    echo "         gawk, samtools, modwt, bedops (bam2bed,bedmap,sort-bed,starch,unstarch) on path."
    exit -1; 
fi
bam_file=$1     # Bam file on which hotspot will be run.
chrom_sizes=$2  # Chrom sizes file that matches the reference genome which bam reads were aligned to.
if [ $# -eq 3 ]; then
    blacklist_bed=$3  # List of regions to exclude from hotspot calling.
    blacklist="-e $blacklist_bed"
else
    blacklist=""
fi

bam_root=${bam_file%.bam}
hotspot_root="${bam_root}_hotspots"
peaks_root="${bam_root}_peaks"
density_root="${bam_root}_density"
### Temporary for debugging
cutcounts_root="${bam_root}_cutcounts"
allcalls_root="${bam_root}_allcalls"
### Temporary for debugging
echo "-- output: '${hotspot_root}.bed', '${peaks_root}.bed', and '${density_root}.bed'"

echo "-- Convert chrom.sizes to bed format..."
cat $chrom_sizes | awk '{printf "%s\t0\t%s\n",$1,$2}' > chrom_sizes.bed

echo "-- Running hotspot2.sh..."
set -x
hotspot2.sh -s 12345 $blacklist -c chrom_sizes.bed $bam_file out/
set +x

if [ -s out/*.hotspots.fdr*.starch ]; then
    echo "-- Converting hotspots to bed and bigBed..."
    set -x
    mv out/*.hotspots.fdr*.starch ${hotspot_root}.starch
    unstarch ${hotspot_root}.starch > ${hotspot_root}.bed
    touch ${hotspot_root}.bb  # First round don't press your luck
    #bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${hotspot_root}.bed $chrom_sizes ${hotspot_root}.bb
    ### currently: 'chr1    10249   10451   id-1    3.14448'
    grep "^chr" ${hotspot_root}.bed | wc -l > ${hotspot_root}_count.txt
    set +x
else
    ### Temporary for debugging
    touch ${hotspot_root}.bed
    touch ${hotspot_root}.bb
    wc -l ${hotspot_root}.bed > ${hotspot_root}_count.txt
fi

if [ -s out/*.peaks.starch ]; then
    echo "-- Converting peaks to bed and bigBed..."
    set -x
    mv out/*.peaks.starch ${peaks_root}.starch
    unstarch ${peaks_root}.starch > ${peaks_root}.bed
    touch ${peaks_root}.bb  # First round don't press your luck
    #bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${peaks_root}.bed $chrom_sizes ${peaks_root}.bb
    grep "^chr" ${peaks_root}.bed | wc -l > ${peaks_root}_count.txt
    set +x
else
    ### Temporary for debugging
    touch ${peaks_root}.bed
    touch ${peaks_root}.bb
    wc -l ${peaks_root}.bed > ${peaks_root}_count.txt
fi

# TODO: Expect bedGraph?
if [ -s out/*.density.starch ]; then
    echo "-- Converting density to bedGrah and bigWig..."
    set -x
    mv out/*.density.starch ${density_root}.starch
    unstarch ${density_root}.starch > ${density_root}.bed
    touch ${density_root}.bw  # First round don't press your luck
    #bedGraphToBigWig ${density_root}.bed $chrom_sizes ${density_root}.bw
    grep "^chr" ${density_root}.bed | wc -l > ${density_root}_count.txt  ### Temporary should be bigWig?
    set +x
else
    ### Temporary for debugging
    touch ${density_root}.bed
    touch ${density_root}.bw
    wc -l ${density_root}.bed > ${density_root}_count.txt
fi

### Temporary for debugging
if [ -s out/*.cutcounts.starch ]; then
    echo "-- Saving cutcounts..."
    set -x
    mv out/*.cutcounts.starch ${cutcounts_root}.starch
    unstarch ${cutcounts_root}.starch | grep "^chr" | wc -l > ${cutcounts_root}_count.txt
    set +x
else
    touch ${cutcounts_root}.starch
    wc -l ${cutcounts_root}.starch > ${cutcounts_root}_count.txt
fi

### Temporary for debugging
if [ -s out/*.allcalls.starch ]; then
    echo "-- Saving allcalls..."
    set -x
    mv out/*.allcalls.starch ${allcalls_root}.starch
    unstarch ${allcalls_root}.starch | grep "^chr" | wc -l > ${allcalls_root}_count.txt
    set +x
else
    touch ${allcalls_root}.starch
    wc -l ${allcalls_root}.starch > ${allcalls_root}_count.txt
fi
echo "-- The results..."
ls -l ${hotspot_root}*
ls -l ${peaks_root}*
ls -l ${density_root}*
ls -l ${cutcounts_root}* ### Temporary for debugging
ls -l ${allcalls_root}* ### Temporary for debugging
df -k .


