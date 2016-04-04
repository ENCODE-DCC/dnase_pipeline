#!/bin/bash -e

if [ $# -lt 6 ] || [ $# -gt 7 ]; then
    echo "usage v1: dnase_hotspot.sh <aligned.bam> <chrom.sizes> <blacklist.bed> <hotspot_root> <peaks_root> <density_root> [<allcalls_root>]"
    echo "Calls peaks hotspot2.  Is independent of DX and encodeD."
    echo "Requires hotspot2 (hotspot2,hotspot2.sh,tallyCountsInSmallWindows,cutcounts.bash,density-peaks.bash,bed_exclude.py),"
    echo "         bedops (bam2bed,bedmap,sort-bed,starch,unstarch), modwt, mawk, samtools, bedToBigBed, bedGraphToBigWig,"
    echo "         and pigz on path."
    exit -1; 
fi
bam_file=$1      # Bam file on which hotspot will be run.
chrom_sizes=$2   # Chrom sizes file that matches the reference genome which bam reads were aligned to.
blacklist_bed=$3 # List of regions to exclude from hotspot calling. (Optional: 'na' for no blacklist).
hotspot_root=$4  # Put hotspot results into ${hotspot_root}.bed.gz, ${hotspot_root}.bb, ${hotspot_root}_count.txt, and ${hotspot_root}_SPOT.txt
peaks_root=$5    # Put peak results into ${peaks_root}.bed.gz, ${peaks_root}.bb, and ${peaks_root}_count.txt
density_root=$6  # Put density results into ${density_root}.bw
echo "-- output: '${hotspot_root}.bed.gz/.bb/_count.txt/_SPOT.txt', '${peaks_root}.bed.gz/.bb/_count.txt', '${density_root}.bw'"
if [ $# -eq 7 ]; then
    allcalls_root=$7 # (optional) Put all_calls results into ${allcalls_root}.bed.gz and ${allcalls_root}_count.txt
    echo "           and '${allcalls_root}.bed.gz/_count.txt'"
fi

echo "-- Convert chrom.sizes to bed format..."
cat $chrom_sizes | awk '{printf "%s\t0\t%s\n",$1,$2}' | sort-bed - > chrom_sizes.bed

blacklist_param=""
if [ -s $blacklist_bed ]; then
    blacklist_param="-e $blacklist_bed"
fi

echo "-- Running hotspot2.sh..."
set -x
hotspot2.sh -s 12345 $blacklist_param -c chrom_sizes.bed $bam_file out/
set +x

echo "-- houtspot2.sh out/..."
ls -l out

echo "-- Converting hotspots to bed.gz and bigBed..."
set -x
mv out/*.hotspots.fdr*.starch ${hotspot_root}.starch
unstarch ${hotspot_root}.starch > ${hotspot_root}.bed
head ${hotspot_root}.bed # FOR DEBUGGING
touch ${hotspot_root}.bb  # First round don't press your luck
#bedToBigBed -as=/usr/bin/broadPeak.as -type=bed6+3 ${hotspot_root}.bed $chrom_sizes ${hotspot_root}.bb
grep "^chr" ${hotspot_root}.bed | wc -l > ${hotspot_root}_count.txt
pigz ${hotspot_root}.bed
set +x

echo "-- Converting narrowpeaks to bed.gz and bigBed..."
set -x
mv out/*.peaks.narrowpeaks.starch ${peaks_root}.starch
unstarch ${peaks_root}.starch > ${peaks_root}.bed
#head ${peaks_root}.bed # FOR DEBUGGING
bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${peaks_root}.bed $chrom_sizes ${peaks_root}.bb
grep "^chr" ${peaks_root}.bed | wc -l > ${peaks_root}_count.txt
pigz ${peaks_root}.bed
set +x

echo "-- Saving density bigWig..."
set -x
mv out/*.density.bw ${density_root}.bw
set +x

echo "-- Saving SPOT..."
set -x
mv out/*.SPOT.txt ${hotspot_root}_SPOT.txt
set +x
echo "-- ----- SPOT score: -----"
cat ${hotspot_root}_SPOT.txt

if [ $# -eq 7 ]; then
    echo "-- Counting allcalls..."
    set -x
    mv out/*.allcalls.starch ${allcalls_root}.starch
    unstarch ${allcalls_root}.starch > ${allcalls_root}.bed
    grep "^chr" ${allcalls_root}.bed | wc -l > ${allcalls_root}_count.txt
    pigz ${allcalls_root}.bed
    set +x
fi

echo "-- The results..."
ls -l ${hotspot_root}*
ls -l ${peaks_root}*
ls -l ${density_root}*
if [ $# -eq 7 ]; then
    ls -l ${allcalls_root}*
fi
df -k .


