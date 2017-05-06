#!/bin/bash -e

if [ $# -lt 3 ] ||  [ $# -gt 4 ]; then
    echo "usage v1: dnase_density.sh <filtered.bam> <chrom.sizes> <density_root> [<chrom_buckets.bed>]"
    echo "Filters single-end aligned reads for DNase.  Is independent of DX and encodeD."
    echo "Requires gawk, samtools, bam2bed, convert2bed, bedmap, sort-bed and wigToBegWig on path."
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

if [ "chrom_buckets_root" != "" ]; then
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
    gawk -v w=${WIN} '{print $$1"\t0\t"$$2-w}' $chrom_sizes | sort-bed - \
      | gawk -v binI=${BINI} -v win=${WIN} '{ for(i = $$2 + win; i < $$3; i += binI) { print $$1"\t"i - win"\t"i + win }}' \
      > $chrom_buckets
    set +x
    ls -l $chrom_buckets
fi

echo "-- Making genome-wide density bed file..."
set -x
bam2bed -d \
	  < ${filtered_bam} \
	  | cut -f1-6 \
    | gawk '{ if( $$6=="+" ){ s=$$2; e=$$2+1 } else { s=$$3-1; e=$$3 } print $$1 "\t" s "\t" e "\tid\t" 1 }' \
    | sort-bed - \
    > ${filtered_bam_root}.bed
set +x
ls -l ${filtered_bam_root}.bed

echo "-- Making bucketted raw (unnormalized) density bed file..."
set -x
cat $chrom_buckets \
    | bedmap --faster --echo --count --delim "\t" - ${filtered_bam_root}.bed \
    | gawk -v binI=${BINI} -v win=${WIN} \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print $$1 "\t" $$2 + shiftFactor "\t" $$3-shiftFactor "\tid\t" i $$4}' \
    > ${raw_density_root}.bed
set +x
ls -l ${raw_density_root}.bed

echo "-- Making normalized density bed file..."
set -x
cat ${raw_density_root}.bed | \
  gawk -v allcounts=`samtools view -c ${filtered_bam}` \
      -v extranuclear_counts=`samtools view -c ${filtered_bam} chrM chrC` \
      -v scale=${SCALE} 'BEGIN{ tagcount=allcounts-extranuclear_counts }{z=$$5; n=(z/tagcount)*scale; print $$1 "\t" $$2 "\t" $$3 "\t" $$4 "\t" n }' \
  | tee ${density_root}.bed \
  | starch - > ${density_root}.starch
set +x
ls -l ${density_root}.*
	  
echo "-- Making wig from bed file..."
set -x
cat ${density_root}.bed \
  | gawk -v s=0 '{ if( $$5 != s ){ print $$0 } }' \
  | gawk 'BEGIN{ OFS = "\t"; chr = ""}{ if( chr == "" ){ chr=$$1; print "variableStep chrom="chr" span=20" } if( $$1 == chr ){ print $$2, $$5 } else { chr=$$1; print "variableStep chrom="chr" span=20"; print $$2, $$5} }' \
  > ${density_root}.wig
set +x
ls -l ${density_root}.*
	  
echo "-- Making BigWig from wig file..."
set -x
wigToBigWig -clip ${density_root}.wig ${chrom_sizes} ${density_root}.bw
set -x

# TODO: any qc ??
#echo "-- Collect bam stats..."
#set -x
#samtools flagstat $unfiltered_bam > ${unfiltered_bam_root}_flagstat.txt
#samtools flagstat ${filtered_bam_root}.bam > ${filtered_bam_root}_flagstat.txt
#samtools stats ${filtered_bam_root}.bam > ${filtered_bam_root}_samstats.txt
#head -3 ${filtered_bam_root}_samstats.txt
#grep ^SN ${filtered_bam_root}_samstats.txt | cut -f 2- > ${filtered_bam_root}_samstats_summary.txt
#set +x

echo "-- The results..."
ls -l ${density_root}*
df -k .

