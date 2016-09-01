#!/bin/bash

usage() {
  cat >&2 <<__EOF__
Usage:  "$0" [options] -c CHROM_SIZES -o OUTFILE

NOTE: This script only needs to be run once per *genome*, not once per *sample*,
      unless a change to the NEIGHBORHOOD_SIZE or MAPPABLE_REGIONS is desired.
      MAPPABLE_REGIONS, when supplied, should not contain any "blacklist" regions
      (problematic satellite repeats, etc.).

Options:
    -h                    Show this helpful help

    -c CHROM_SIZES        BED or starch file of chromosome sizes, with column 2 set to zeroes. Mandatory.
    -o OUTFILE            Output file name. Mandatory. If it doesn't end in .starch, .starch will be appended.

    -M MAPPABLE_REGIONS   BED or starch file of mappable regions, with "blacklist" subtracted when appropriate.
    -n NEIGHBORHOOD_SIZE  Local neighborhood radius in bp (default = 100, yielding a 201-bp window)
__EOF__
  exit 2
}

# Input parsing
CHROM_SIZES=
MAPPABLE_REGIONS=
OUTFILE=
HALF_WINDOW_SIZE=100

while getopts 'hc:M:o:n:' opt; do
  case "$opt" in
    h) usage ;;
    c) CHROM_SIZES=$OPTARG ;;
    M) MAPPABLE_REGIONS=$OPTARG ;;
    o) OUTFILE=$OPTARG ;;
    n) HALF_WINDOW_SIZE=$OPTARG ;;
  esac
done

if [ "$CHROM_SIZES" == "" ]; then
  echo -e "Error:  Required argument -c CHROM_SIZES was not provided."
  usage
fi

if [ ! -s "$CHROM_SIZES" ]; then
  echo -e "Error:  CHROM_SIZES file \"$CHROM_SIZES\" was not found, or is empty."
  usage
fi

if [ "$MAPPABLE_REGIONS" != "" ] && [ ! -s "$MAPPABLE_REGIONS" ]; then
  echo -e "Error:  MAPPABLE_REGIONS file \"$MAPPABLE_REGIONS\" was not found, or is empty."
  usage
fi

if [ "$OUTFILE" == "" ]; then
  echo -e "Error:  Required argument -o OUTFILE was not provided."
  usage
fi

# Force output file name to end in .starch.
if [[ "$OUTFILE" != *.starch ]]; then
  OUTFILE=$OUTFILE.starch
fi

HALFWINSIZE_IS_POSITIVE_INTEGER=$(echo "$HALF_WINDOW_SIZE" | awk '{len=split($0,x,".");if(1==len && int($0)>0){print 1}else{print 0}}')
if [ "$HALFWINSIZE_IS_POSITIVE_INTEGER" != "1" ]; then
  echo -e "Invalid NEIGHBORHOOD_SIZE \"$HALF_WINDOW_SIZE\"; must be a positive integer (default = 100)."
  usage
fi

# Get all sites (1bp each) that can be viable centers of windows in which we'll want to tally cut counts.
# This means all mappable sites that are not within a half-window's width
# of any unmappable region whose width is >= the half-window width.
if [ "$MAPPABLE_REGIONS" != "" ]; then
  bedops --range "${HALF_WINDOW_SIZE}:-${HALF_WINDOW_SIZE}" -u "$CHROM_SIZES" \
    | bedops -d - "$MAPPABLE_REGIONS" \
    | awk -v "t=$HALF_WINDOW_SIZE" '{if($3-$2>=t){beg=$2-t;if(beg<0){beg=0}print $1"\t"beg"\t"$3+t}}' \
    | bedops -d "$MAPPABLE_REGIONS" - \
    | bedops -w - \
    | starch - \
      >"$OUTFILE"
else
  bedops --range "${HALF_WINDOW_SIZE}:-${HALF_WINDOW_SIZE}" -u "$CHROM_SIZES" \
    | bedops -w - \
    | starch - \
      >"$OUTFILE"
fi
