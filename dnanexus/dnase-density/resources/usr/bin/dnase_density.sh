#!/bin/bash -e

if [ $# -lt 3 ] ||  [ $# -gt 4 ]; then
    echo "usage v1: dnase_density.sh <filtered.bam> <chrom.sizes> <density_root> [<chrom_buckets.bed>]"
    echo "Filters single-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires gawk, samtools, bam2bed, convert2bed, bedmap, sort-bed and bedGraphToBigWig on path."
    exit -1; 
fi
filtered_bam=$1  # filtered bam file.
chrom_sizes=$2   # chrom sizes file that matches the reference genome which bam reads were aligned to.
density_root=$3  # root name for density file output (e.g. "out" will create "out.bw") 
if [ $# -eq 4 ]; then
    chrom_buckets=$4 # chrom_buckets file in either bed.gz or starch format
else
    chrom_buckets=""
fi

filtered_bam_root=${filtered_bam%.bam}
raw_density_root=${filtered_bam_root}_raw
echo "-- Normalized density file will be: '${density_root}.bw'"

BINI=20
WIN=75
SCALE=1000000

if [ "$chrom_buckets" != "" ]; then
    chrom_buckets_root=${chrom_buckets%.bed.gz}
    if [ "$chrom_buckets_root" != "$chrom_buckets" ]; then
        echo "-- gunzipping ${chrom_buckets} file..."
        set -x
        gunzip $chrom_buckets
        set +x
        chrom_buckets=${chrom_buckets_root}.bed
    else
        chrom_buckets_root=${chrom_buckets%.starch}
        if [ "$chrom_buckets_root" != "$chrom_buckets" ]; then
            echo "-- unstarching ${chrom_buckets} file..."
            set -x
            unstarch $chrom_buckets > ${chrom_buckets_root}.bed
            set +x
            chrom_buckets=${chrom_buckets_root}.bed
        fi
    fi
else
    chrom_buckets=chrom_buckets.bed
    echo "-- Making ${chrom_buckets} file..."
    set -x
    gawk -v win=${WIN} '{print $1"\t0\t"$2-win}' $chrom_sizes | sort-bed - \
      | gawk -v binI=${BINI} -v win=${WIN} '{ for(i = $2 + win; i < $3; i += binI) { print $1"\t"i - win"\t"i + win }}' \
      > $chrom_buckets
    set +x
fi
ls -l $chrom_buckets

echo "-- Making genome-wide density bed file..."
set -x
bam2bed -d \
	  < ${filtered_bam} \
	  | cut -f1-6 \
    | gawk '{ if( $6=="+" ){ s=$2; e=$2+1 } else { s=$3-1; e=$3 } print $1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > ${filtered_bam_root}.bed
set +x
ls -l ${filtered_bam_root}.bed

echo "-- Making windowed raw (unnormalized) density bed file..."
set -x
cat $chrom_buckets \
    | bedmap --faster --echo --count --delim "\t" - ${filtered_bam_root}.bed \
    | gawk -v binI=${BINI} -v win=${WIN} \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print $1 "\t" $2 + shiftFactor "\t" $3-shiftFactor "\tid\t" i $4}' \
    > ${raw_density_root}.bed
set +x
ls -l ${raw_density_root}.bed

echo "-- Making normalized density bed file..."
set -x
cat ${raw_density_root}.bed | \
  gawk -v allcounts=`samtools view -c ${filtered_bam}` \
      -v extranuclear_counts=`samtools view -c ${filtered_bam} chrM chrC` \
      -v scale=${SCALE} 'BEGIN{ tagcount=allcounts-extranuclear_counts }{z=$5; n=(z/tagcount)*scale; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" n }' \
  > ${density_root}.bed
starch ${density_root}.bed > ${density_root}.starch
set +x
ls -l ${density_root}.*
	  
echo "-- Making BigWig from bed file..."
set -x
cat ${density_root}.bed | cut -f1,2,3,5 > density.bedGraph
bedGraphToBigWig density.bedGraph $chrom_sizes ${density_root}.bw
set +x
ls -l ${density_root}.bw
	  
# qc anyone?  Not much value.
echo "-- Count density intervals..."
set -x
grep "^chr" ${density_root}.bed | wc -l > ${density_root}_count.txt
set +x

echo "-- The results..."
ls -l ${density_root}*
df -k .

