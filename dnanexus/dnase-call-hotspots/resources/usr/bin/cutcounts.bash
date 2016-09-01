#!/bin/bash
# This is a modified version of the cutcounts script used in stampipes
set -x -e -o pipefail

if [[ $# != 5 && $# != 6 ]]; then
  echo "Usage: $0 in.bam cutcounts.starch fragments.starch totalcuts.txt chromSizes [mappableRegions.starch]" >&2
  echo -e "where the first three arguments after the BAM filename contain names of your choosing"
  echo -e "for output files that will be created, and chromSizes is a .bed or .starch file of chromosome sizes,"
  echo -e "with the start of each chromosome (column 2) set to 0."
  echo -e "A file of mappable regions can optionally be supplied as a 4th argument; this is recommended."
  exit 2
fi

bam=$1
CUTCOUNTS=$2
FRAGMENTS=$3
TOTALCUTSFILE=$4
CHROM_SIZES=$5
MAPPABLE_REGIONS=""
if [[ $# == 6 ]]; then
  MAPPABLE_REGIONS=$6
  if [ ! -s "$MAPPABLE_REGIONS" ]; then
    echo -e "Error:  File \"$MAPPABLE_REGIONS\" supplied to $0 was not found, or is empty."
    exit 2
  fi
fi
if [ ! -s "$bam" ]; then
  echo -e "Error:  Alignment file \"$bam\" supplied to $0 was not found, or is empty."
  exit 2
fi
if [ ! -s "$CHROM_SIZES" ]; then
  echo -e "Error:  Chromosome sizes file \"$CHROM_SIZES\" supplied to $0 was not found, or is empty."
  exit 2
fi

clean=0
if [[ -z "$TMPDIR" ]]; then
  TMPDIR=$(mktemp -d)
  clean=1
fi

# Prefer mawk, if installed
AWK_EXE=$(which mawk 2>/dev/null || which awk)

# temp files
FRAGMENTSTMP="$TMPDIR/fragments.bed"

# Create cut counts and fragments if they don't exist
if [[ ! -s "$CUTCOUNTS" ]]; then

  if [ "$MAPPABLE_REGIONS" == "" ]; then
    time bam2bed --do-not-sort <"$bam" \
      | "$AWK_EXE" -v "fragmentfile=$FRAGMENTSTMP"
    'BEGIN { FS="\t"; OFS=FS }
          {
            strand=$6;
            read_start=$2;
            read_end=$3;
            read_id=$1;
            flag=$7;
            tlen=$11;
            if( strand == "+" ) {
              cut_start = read_start;
              cut_end = read_start + 1;
            } else {
              cut_start= read_end;
              cut_end = read_end + 1;
            }
            print read_id, cut_start, cut_end;
            if (tlen > 0) {
              fragment_end = read_start + tlen;
              print read_id, read_start, fragment_end > fragmentfile;
            }
          }' \
      | sort-bed --max-mem 8G - \
      | uniq -c \
      | "$AWK_EXE" '{ print $2"\t"$3"\t"$4"\tid-"NR"\t"$1 }' \
      | bedops -e 1 - "$CHROM_SIZES" \
      | "$AWK_EXE" -v "totalcutsfile=$TOTALCUTSFILE" 'BEGIN{sum=0}{sum += $5; print $0}END{print sum > totalcutsfile}' \
      | starch - \
        >"$CUTCOUNTS"

  else

    time bam2bed --do-not-sort <"$bam" \
      | "$AWK_EXE" -v "fragmentfile=$FRAGMENTSTMP" \
        'BEGIN { FS="\t"; OFS=FS }
          {
            strand=$6;
            read_start=$2;
            read_end=$3;
            read_id=$1;
            flag=$7;
            tlen=$11;
            if( strand == "+" ) {
              cut_start = read_start;
              cut_end = read_start + 1;
            } else {
              cut_start= read_end;
              cut_end = read_end + 1;
            }
            print read_id, cut_start, cut_end;
            if (tlen > 0) {
              fragment_end = read_start + tlen;
              print read_id, read_start, fragment_end > fragmentfile;
            }
          }' \
      | sort-bed --max-mem 8G - \
      | uniq -c \
      | "$AWK_EXE" '{ print $2"\t"$3"\t"$4"\tid-"NR"\t"$1 }' \
      | bedops -e 1 - "$CHROM_SIZES" \
      | bedops -e -1 - "$MAPPABLE_REGIONS" \
      | "$AWK_EXE" -v "totalcutsfile=$TOTALCUTSFILE" 'BEGIN{sum=0}{sum += $5; print $0}END{print sum > totalcutsfile}' \
      | starch - \
        >"$CUTCOUNTS"

  fi

  if [[ ! -e "$FRAGMENTSTMP" ]]; then
    touch "$FRAGMENTSTMP"
  fi

  if [ "$MAPPABLE_REGIONS" == "" ]; then
    sort-bed --max-mem 8G "$FRAGMENTSTMP" \
      | bedops -e 1 - "$CHROM_SIZES" \
      | starch - >"$FRAGMENTS"
  else
    sort-bed --max-mem 8G "$FRAGMENTSTMP" \
      | bedops -e 1 - "$CHROM_SIZES" \
      | bedops -e 1 - "$MAPPABLE_REGIONS" \
      | starch - >"$FRAGMENTS"
  fi

  rm -f "$FRAGMENTSTMP"

fi

if [ $clean != 0 ]; then
  rm -rf "$TMPDIR"
fi

exit 0
