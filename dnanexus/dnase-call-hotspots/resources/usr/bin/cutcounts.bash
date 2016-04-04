#!/bin/bash
# This is a modified version of the cutcounts script used in stampipes
set -x -e -o pipefail

if [[ $# != 4 ]] ; then
  echo "Usage: $0 in.bam cutcounts.starch fragments.starch" >&2
  exit 2
fi

bam=$1
CUTCOUNTS=$2
FRAGMENTS=$3
TOTALCUTSFILE=$4

name=$(basename $bam .bam)
outputdir="$(dirname $CUTCOUNTS)"

clean=0
if [[ -z "$TMPDIR" ]] ;then
  TMPDIR=$(mktemp -d)
  clean=1
fi

# Prefer mawk, if installed
AWK_EXE=$(which mawk 2>/dev/null || which awk)

# temp files
FRAGMENTSTMP="$TMPDIR/fragments.bed"

# Create cut counts and fragments if they don't exist
if [[  ! -s "$CUTCOUNTS" ]]; then

  time bam2bed --do-not-sort < "$bam" \
    | "$AWK_EXE" -v fragmentfile=$FRAGMENTSTMP -v totalcutsfile=$TOTALCUTSFILE \
        'BEGIN { FS="\t"; OFS=FS } ; { \
          strand=$6; \
          read_start=$2; \
          read_end=$3; \
          read_id=$1; \
          flag=$7; \
          tlen=$11; \
          if( strand == "+" ) { \
            cut_start = read_start; \
            cut_end = read_start + 1; \
          } else { \
            cut_start= read_end; \
            cut_end = read_end + 1; \
          } \
          print read_id, cut_start, cut_end; \
          if (tlen > 0) { \
            fragment_end = read_start + tlen; \
            print read_id, read_start, fragment_end > fragmentfile; \
          } \
        } END { print NR > totalcutsfile }' \
    | sort-bed --max-mem 8G - \
    | uniq -c \
    | "$AWK_EXE" '{ print $2"\t"$3"\t"$4"\tid-"NR"\t"$1 }' \
    | starch - \
    > "$CUTCOUNTS"

  if [[  ! -e "$FRAGMENTSTMP" ]]; then
    touch "$FRAGMENTSTMP"
  fi
  sort-bed --max-mem 8G "$FRAGMENTSTMP" | starch - > "$FRAGMENTS"
  rm -f "$FRAGMENTSTMP"

fi

if [ $clean != 0 ] ; then
  rm -rf $TMPDIR
fi

exit 0
