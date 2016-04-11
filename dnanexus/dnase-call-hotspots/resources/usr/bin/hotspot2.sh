#! /bin/bash
set -x -e -o pipefail

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
    -m MIN_HOTSPOT_WIDTH  Minimum hotspot width allowed (50)
    -p PVAL_FDR           Number of p-values to use for FDR (1000000)
    -f HOTSPOT_THRESHOLD  False-discovery rate to use for hotspot filtering (0.05)
    -F SITECALL_THRESHOLD False-discovery rate to use for site-call filtering (0.10)
    -O                    Use non-overlapping windows (advanced option)

    -s SEED               Set this to an integer for repeatable results

    Both the exclude file and chromosome sizes file should be in bed or starch
    format.

    Neighborhood and window sizes are specified as the distance from the edge
    to the center - i.e, a 100bp neighborhood size is a 201bp window.

    The site-call False Discovery Rate threshold (-F) must be greater than or
    equal to the hotspot FDR threshold.

    Using non-overlapping windows is not recommended for most users.

__EOF__
    exit 2
}


EXCLUDE_THESE_REGIONS="/dev/null"
CHROM_SIZES=""
SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=100 # i.e., 201bp regions
BACKGROUND_WINDOW_SIZE=50001 # i.e., +/-25kb around each position
MIN_HOTSPOT_WIDTH=50
PVAL_DISTN_SIZE=1000000
OVERLAPPING_OR_NOT="overlapping"
HOTSPOT_FDR_THRESHOLD="0.05"
CALL_THRESHOLD="$HOTSPOT_FDR_THRESHOLD"
SEED=$RANDOM

while getopts 'hc:e:f:F:m:n:p:s:w:O' opt ; do
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
      HOTSPOT_FDR_THRESHOLD=$OPTARG
      ;;
    F)
      CALL_THRESHOLD=$OPTARG
      ;;
    m)
      MIN_HOTSPOT_WIDTH=$OPTARG
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

# Check to make sure Hotspot FDR <= site-call FDR
if awk '{exit $1>$2?0:1}' <<< "$HOTSPOT_FDR_THRESHOLD $CALL_THRESHOLD"; then
  echo "Hotspot FDR threshold (-f) cannot be greater than site-calling threshold (-F)" >&2
  exit 2
fi

BAM=$1
OUTDIR=$2

echo "checking system for modwt, BEDOPS, samtools..."
if ! which modwt &>/dev/null; then
  echo "Could not find modwt!"
  exit -1
elif ! which bedGraphToBigWig &>/dev/null; then
  echo "Could not find bedGraphToBigWig!"
  exit -1
elif ! which bedmap &>/dev/null; then
  echo "Could not find BEDOPS (bedmap)!"
  exit -1
elif ! which samtools &>/dev/null; then
  echo "Could not find samtools!"
  exit -1
elif ! which hotspot2 &>/dev/null; then
  echo "Could not find hotspot2!"
  exit -1
elif ! which tallyCountsInSmallWindows &>/dev/null; then
  echo "Could not find tallyCountsInSmallWindows!"
  exit -1
fi
WAVELETS_EXE=$(which modwt)
CUTCOUNT_EXE="$(dirname "$0")/cutcounts.bash"
DENSPK_EXE="$(dirname "$0")/density-peaks.bash"
EXCLUDE_EXE="$(dirname "$0")/bed_exclude.py"
COUNTING_EXE=tallyCountsInSmallWindows
HOTSPOT_EXE=hotspot2

# Prefer mawk, if installed
AWK_EXE=$(which mawk 2>/dev/null || which awk)


mkdir -p "$OUTDIR"

base="$OUTDIR/$(basename "$BAM" .bam)"

HOTSPOT_OUTFILE="$base.hotspots.fdr$HOTSPOT_FDR_THRESHOLD.starch"
CUTCOUNTS="$base.cutcounts.starch"
FRAGMENTS_OUTFILE="$base.fragments.sorted.starch"
TOTALCUTS_OUTFILE="$base.cleavage.total"
OUTFILE="$base.allcalls.starch"
DENSITY_OUTFILE="$base.density.starch"
DENSITY_BW="$base.density.bw"
PEAKS_OUTFILE="$base.peaks.starch"
SPOT_SCORE_OUTFILE="$base.SPOT.txt"


clean=0
if [[ -z "$TMPDIR" ]] ;then
  TMPDIR=$(mktemp -d)
  clean=1
fi

echo "Cutting..."
bash "$CUTCOUNT_EXE" "$BAM" "$CUTCOUNTS" "$FRAGMENTS_OUTFILE" "$TOTALCUTS_OUTFILE"

# don't unstarch $CUTCOUNTS and feed to $COUNTING_EXE since things like chrM may be in $CUTCOUNTS but not $CHROM_SIZES
# run $CHROM_SIZES through bedops --ec -u for error checking
echo "Running hotspot2..."
bedops --ec -u "$CHROM_SIZES" \
    | bedops --ec -e 1 "$CUTCOUNTS" - \
    | "$COUNTING_EXE" "$SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE" "$OVERLAPPING_OR_NOT" "reportEachUnit" "$CHROM_SIZES" \
    | python "$EXCLUDE_EXE" "$EXCLUDE_THESE_REGIONS" \
    | "$HOTSPOT_EXE" --fdr_threshold "$CALL_THRESHOLD" --background_size="$BACKGROUND_WINDOW_SIZE" --num_pvals="$PVAL_DISTN_SIZE" --seed="$SEED" \
    | starch - \
    > "$OUTFILE"

# We report the largest -log10(P) observed at any bp of a hotspot
# as the "score" of that hotspot, where P is the site-specific P-value.
# P-values of 0 will be encountered, and we don't want to do log(0).
# Nonzero P-values as low as 1e-308 have been seen during testing.
# We choose to cap all P-values at 1e-100, i.e. -log10(P) = 100.
# The constant c below converts from natural logarithm to log10.

echo "Thresholding..."
unstarch "$OUTFILE" \
    | "$AWK_EXE" -v "threshold=$HOTSPOT_FDR_THRESHOLD" '($6 <= threshold)' \
    | bedops -m - \
    | "$AWK_EXE" -v minW=$MIN_HOTSPOT_WIDTH 'BEGIN{FS="\t";OFS="\t"}{ \
          chrR=$1;begPosR=$2;endPosR=$3;widthR=endPosR-begPosR; \
          if(NR>1) { \
            if (chrR == chrL) { \
              distLR = begPosR - endPosL; \
            } \
            else { \
              distLR=999999999; \
            } \
            if (distLR > minW) { \
              if (widthL >= minW) { \
                print chrL, begPosL, endPosL; \
              } \
            } \
            else { \
              if (widthL < minW || widthR < minW) { \
                begPosR = begPosL; \
                widthR = endPosR - begPosR; \
              } \
              else { \
                print chrL, begPosL, endPosL; \
              } \
            } \
          } \
          chrL = chrR; \
          begPosL = begPosR; \
          endPosL = endPosR; \
          widthL = widthR; \
        } END{if(widthL >= minW){print chrL, begPosL, endPosL}}' \
    | sort-bed --max-mem 1G - \
    | bedmap --faster --sweep-all --delim "\t" --echo --min - "$OUTFILE" \
    | bedmap --faster --sweep-all --delim "\t" --echo --count - "$CUTCOUNTS" \
    | "$AWK_EXE" 'BEGIN{OFS="\t";c=-0.4342944819}
        {
          if($4>1e-100) {
            print $1, $2, $3, "id-"NR, $5, ".","-1","-1", c*log($4)
          } else {
            print $1, $2, $3, "id-"NR, $5, ".","-1","-1", "100"
          }
        }' \
     | starch - \
     > "$HOTSPOT_OUTFILE"

unstarch "$HOTSPOT_OUTFILE" \
  | awk 'BEGIN {OFS="\t"} { if ( $5 > 1000 ) { $5 = 1000 } print; }' \
  | starch - \
  > ${HOTSPOT_OUTFILE/.starch/.broadpeaks.starch}

echo "Calculating SPOT score..."
num_cleaves=$(cat "$TOTALCUTS_OUTFILE")
cleaves_in_hotspots=$(bedops --ec -e -1 "$CUTCOUNTS" "$HOTSPOT_OUTFILE" | awk 'BEGIN{s=0} {s+=$5} END {print s}')
echo "scale=4; $cleaves_in_hotspots / $num_cleaves" \
  | bc \
  > "$SPOT_SCORE_OUTFILE"

echo "Generating peaks and density..."
bash "$DENSPK_EXE" "$TMPDIR" "$WAVELETS_EXE" "$CUTCOUNTS" "$HOTSPOT_OUTFILE" "$CHROM_SIZES" "$DENSITY_OUTFILE" "$PEAKS_OUTFILE"

TMPFRAGS="$(mktemp -t fragsXXXXX)"
unstarch "$DENSITY_OUTFILE" | cut -f1,2,3,5 > "$TMPFRAGS"
bedGraphToBigWig \
  "$TMPFRAGS" \
  <(cut -f1,3 "$CHROM_SIZES") \
  "$DENSITY_BW"

echo "Done!"

rm -f "$TMPFRAGS"
if [[ $clean != 0 ]] ; then
  rm -rf "$TMPDIR"
fi

exit 0
