#! /bin/bash

usage(){
  cat >&2 <<__EOF__
Usage:  "$0" [options] in.bam outdir

Options:
    -h                    Show this helpful help

  Mandatory options:
    -c CHROM_SIZES_FILE   The length of each chromosome, in BED or starch format
                          Must include all chromosomes present in BAM input

  Recommended options:
    -e EXCLUDE_FILE       Exclude these regions from analysis

  Optional options:
    -n NEIGHBORHOOD_SIZE  Local neighborhood size (100)
    -w WINDOW_SIZE        Background region size  (25000)
    -p PVAL_FDR           Number of p-values to use for FDR (1000000)
    -f FDR_THRESHOLD      The false-discovery rate to use for filtering (0.05)
    -O                    Use non-overlapping windows (advanced option)

    -s SEED               Set this to an integer for repeatable results

    Both the exclude file and chromosome sizes file should be in bed or starch
    format.

    Neighborhood and window sizes are specified as the distance from the edge
    to the center - i.e, a 100bp neighborhood size is a 201bp window.

    Using non-overlapping windows is not recommended for most users.

__EOF__
    exit 2
}


EXCLUDE_THESE_REGIONS="/dev/null"
CHROM_SIZES=""
SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=100 # i.e., 201bp regions
BACKGROUND_WINDOW_SIZE=50001 # i.e., +/-25kb around each position
PVAL_DISTN_SIZE=1000000
OVERLAPPING_OR_NOT="overlapping"
FDR_THRESHOLD="0.05"
SEED=""

while getopts 'hc:e:f:m:n:p:s:w:O' opt ; do
  case "$opt" in
    h)
      usage
      ;;
    c)
      CHROM_SIZES=$OPTARG
      ;;
    e)
      EXCLUDE_THESE_REGIONS=$OPTARG
      ;;
    f)
      FDR_THRESHOLD=$OPTARG
      ;;
    n)
      SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=$OPTARG
      ;;
    O)
      OVERLAPPING_OR_NOT="nonoverlapping"
      ;;
    p)
      PVAL_DISTN_SIZE=$OPTARG
      ;;
    s)
      SEED=$OPTARG
      ;;
    w)
      BACKGROUND_WINDOW_SIZE=$(( 2 * OPTARG + 1 ))
      ;;

  esac
done
shift $((OPTIND-1))

if [[ -z "$1" || -z "$2" || -z "$CHROM_SIZES" ]]; then
  usage
fi

BAM=$1
OUTDIR=$2

echo "checking system for modwt, BEDOPS, samtools, ..."
if [ $(which modwt &>/dev/null || echo "$?") ] ; then
  echo "Could not find modwt!"
  exit -1
elif [ $(which bedmap &>/dev/null || echo "$?") ] ; then
  echo "Did not find BEDOPS (bedmap)!"
  exit -1
elif [ $(which samtools &>/dev/null || echo "$?") ] ; then
  echo "Did not find samtools!"
  exit -1
elif [ $(which hotspot2 &>/dev/null || echo "$?") ] ; then
  echo "Did not find hotspot2!"
  exit -1
elif [ $(which tallyCountsInSmallWindows &>/dev/null || echo "$?") ] ; then
  echo "Did not find tallyCountsInSmallWindows!"
  exit -1
fi
WAVELETS_EXE=$(which modwt)
CUTCOUNT_EXE="$(dirname "$0")/cutcounts.bash"
DENSPK_EXE="$(dirname "$0")/density-peaks.bash"
COUNTING_EXE=tallyCountsInSmallWindows
HOTSPOT_EXE=hotspot2

mkdir -p $OUTDIR

HOTSPOT_OUTFILE="$OUTDIR/$(basename "$BAM" .bam).hotspots.fdr"$FDR_THRESHOLD".starch"
CUTCOUNTS="$OUTDIR/$(basename "$BAM" .bam).cutcounts.starch"
OUTFILE="$OUTDIR/$(basename "$BAM" .bam).allcalls.starch"
DENSITY_OUTFILE="$OUTDIR/$(basename "$BAM" .bam).density.starch"
PEAKS_OUTFILE="$OUTDIR/$(basename "$BAM" .bam).peaks.starch"

TMPDIR=${TMPDIR:-$(mktemp -d)}

echo "Cutting..."
bash "$CUTCOUNT_EXE" "$BAM" "$CUTCOUNTS"

# don't unstarch $CUTCOUNTS and feed to $COUNTING_EXE since things like chrM may be in $CUTCOUNTS but not $CHROM_SIZES
echo "Running hotspot2..."
bedops -e 1 "$CUTCOUNTS" "$CHROM_SIZES" \
    | "$COUNTING_EXE" "$SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE" "$OVERLAPPING_OR_NOT" "reportEachUnit" "$CHROM_SIZES" \
    | bedops -n 1 - "$EXCLUDE_THESE_REGIONS" \
    | "$HOTSPOT_EXE" "$BACKGROUND_WINDOW_SIZE" "$PVAL_DISTN_SIZE" $SEED \
    | starch - \
    > "$OUTFILE"


# P-values of 0 will exist, and we don't want to do log(0).
# Roughly 1e-308 is the smallest nonzero P usually seen,
# so we can cap everything at that, or use a different tiny value.
# The constant c below converts from natural logarithm to log10.

echo "Thresholding..."
unstarch "$OUTFILE" \
    | awk -v "threshold=$FDR_THRESHOLD" '($6 <= threshold)' \
    | bedops -m - \
    | bedmap --faster --sweep-all --delim "\t" --echo --min - "$OUTFILE" \
    | awk 'BEGIN{OFS="\t";c=-0.4342944819}
      {
        if($4>1e-308) {
          print $1, $2, $3, "id-"NR, c*log($4)
        } else {
          print $1, $2, $3, "id-"NR, "308"
        }
      }' \
     | starch - \
     > "$HOTSPOT_OUTFILE"

bash "$DENSPK_EXE" "$TMPDIR" "$WAVELETS_EXE" "$CUTCOUNTS" "$HOTSPOT_OUTFILE" "$CHROM_SIZES" "$DENSITY_OUTFILE" "$PEAKS_OUTFILE"

echo "Done!"

exit 0
