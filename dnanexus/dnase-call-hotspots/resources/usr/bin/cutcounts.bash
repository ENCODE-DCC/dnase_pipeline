#!/bin/bash
# This is a modified version of the cutcounts script used in stampipes
set -x -e -o pipefail

if [[ $# != 2 ]] ; then
  echo "Usage: $0 in.bam output.starch" >&2
  exit 2
fi

bam=$1
CUTCOUNTS=$2

name=$(basename $bam .bam)
outputdir="$(dirname $CUTCOUNTS)"

FRAGMENTS=$outputdir/$name.fragments.sorted.starch

clean=0
if [[ -z "$TMPDIR" ]] ;then
  TMPDIR=$(mktemp -d)
  clean=1
fi

# temp files
FRAGMENTSTMP="$TMPDIR/fragments.bed"

# Create cut counts and fragments if they don't exist
if [[  ! -s "$CUTCOUNTS" || ! -s "$FRAGMENTS" ]]; then

  time bam2bed --do-not-sort < "$bam" \
    | awk -v fragmentfile=$FRAGMENTSTMP \
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
        }' \
    | sort-bed --max-mem 8G - \
    | uniq -c \
    | awk '{ print $2"\t"$3"\t"$4"\tid-"NR"\t"$1 }' \
    | starch - \
    > "$CUTCOUNTS"

  ### DCC removed: ###sort-bed --max-mem 8G "$FRAGMENTSTMP" | starch - > "$FRAGMENTS"
  ### DCC added
  if [ -s "$FRAGMENTSTMP" ]; then
    sort-bed --max-mem 8G "$FRAGMENTSTMP" | starch - > "$FRAGMENTS"
  else
    echo "*** $FRAGMENTSTMP is empty!"
  fi
  ### DCC added

  rm -f "$FRAGMENTSTMP"

fi

if [ $clean != 0 ] ; then
  rm -rf $TMPDIR
fi

exit 0
