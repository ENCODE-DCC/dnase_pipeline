#!/bin/bash
set -x -e -o pipefail

log() {
  echo -e "$(date '+%Y-%m-%d %H:%M:%S')\t$*"
}

if [[ $# != 7 ]]; then
  echo "Usage: $0 tmpdir wavelets-binary <tags.starch> <hotspots.starch> <chrom-sizes.starch> <density-out.starch> <peaks-out.starch>" >&2
  exit 2
fi

tmpdir=$1
wavelets=$2
tags=$3
hotspots=$4
chrfile=$5
density=$6
pk=$7

if [ ! -d "$tmpdir" ]; then
  mkdir -p "$tmpdir"
fi

## density params
bins=150
step=20
halfbin=$((bins / 2))
rangepad=$((bins / 2 - step / 2))

## wavelet peakfinding params
waveletlvl=3
filter_type=Haar
boundary_type=reflected

# Prefer mawk, if installed
AWK_EXE=$(which mawk 2>/dev/null || which awk)

log "Calculating densities and peak-finding..."
pkouts=""
densouts=""
for chr in $(awk '{print $1}' "$chrfile"); do
  log "\tProcessing $chr"

  ## Tag density, 150bp window, sliding every 20bp, used for peak-finding and display
  ##  --sweep-all used to prevent a possible broken pipe
  ## Use awk, not sed to change NAN to 0 since a 'chromosome' name could include NAN and then so would change fields 1 and 4
  ## sed could be used with regular expressions, but awk is easy enough here, I think.
  bedops --ec --chop 20 --stagger 20 --chrom "$chr" "$chrfile" \
    | bedmap --faster --sweep-all --chrom "$chr" --range "$rangepad" --delim "\t" --prec 0 --echo --echo-ref-row-id --sum - "$tags" \
    | awk 'BEGIN {OFS="\t"} ; { if ( $(NF) == "NAN" ) { $(NF)=0; } print; }' \
    | starch - \
      >"$tmpdir/.dens.$chr.starch"

  densouts="$densouts $tmpdir/.dens.$chr.starch"

  unstarch "$tmpdir/.dens.$chr.starch" \
    | cut -f5 \
    | $wavelets --level $waveletlvl --to-stdout --boundary $boundary_type --filter $filter_type - \
      >"$tmpdir/.waves"

  ## changed from Bob's which printed out the wavelet smoothed value for a peak, instead
  ## print the density value so things match up with the forced peak-per-hotspots
  ##   bedops -n uses --ec to prevent possible broken pipe
  unstarch "$tmpdir/.dens.$chr.starch" \
    | paste - "$tmpdir/.waves" \
    | "$AWK_EXE" 'BEGIN{incr=0} ; {
            if ( NR > 0 ) {
              if ( incr == 1 ) {
                if ( $6-lastv < 0 ) {
                  print lastl; incr=0;
                }
              } else {
                if ( $6-lastv > 0 ) incr=1;
              }
            }
            lastv=$6; lastl=$0;
          }' \
    | cut -f1-5 \
    | tee "$tmpdir/.wave-pks.$chr" \
    | bedops --ec --chrom "$chr" -n 1 "$hotspots" - \
      >"$tmpdir/.hots-no-pks.$chr"

  ## force a peak call in hotspots with no current peak calls ('peak-per-hotspot')
  ##   --sweep-all to prevent a possible broken pipe; --prec 1 shows decimals for peaks forced by hotspot calls
  bedmap --faster --sweep-all --prec 1 --max-element "$tmpdir/.hots-no-pks.$chr" "$tmpdir/.dens.$chr.starch" \
    | sort-bed - \
    | bedops -u - "$tmpdir/.wave-pks.$chr" \
      >"$tmpdir/.full.pks.$chr"

  pkouts="$pkouts $tmpdir/.full.pks.$chr"

  rm -f "$tmpdir/.waves" "$tmpdir/.wave-pks.$chr" "$tmpdir/.hots-no-pks.$chr"
done

log "Finalizing peaks..."
cat $pkouts \
  | "$AWK_EXE" -v "h=$halfbin" -v "c=$chrfile" \
    'BEGIN {
      OFS="\t";
      while ( (getline line < c) > 0 ) {
        split(line, a, "\t");
        chrom_ends[a[1]] = a[3];
      }
    } ; {
      m=($2+$3)/2; left=m-h; if(left < 0) left=0;
      right=m+h; if(chrom_ends[$1] < right) right=chrom_ends[$1];
      print $1, left, right, $4, $5;
    }' \
    - \
  | bedmap --echo --skip-unmapped --sweep-all --fraction-either 0.25 - "$hotspots" \
  | starch - \
    >"$pk"

unstarch "$pk" \
  | awk -v "h=$halfbin" 'BEGIN {OFS="\t"} ; { print $1, $2, $3, ".", "0", ".", $5, "-1", "-1", h }' \
  | starch - \
    >"${pk/.starch/.narrowpeaks.starch}"

log "Finalizing density..."
starchcat $densouts \
  >"$density"

exit 0
