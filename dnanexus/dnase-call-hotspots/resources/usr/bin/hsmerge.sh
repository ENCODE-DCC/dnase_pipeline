#!/bin/bash

usage(){
  cat >&2 <<__EOF__
Usage:  "$0" [options] input.allcalls.starch hotspots.out.starch

Options:
    -h                  Show this helpful help

    -c CUTCOUNTS        Cut counts file. Mandatory.

    -f FDR_THRESHOLD    Sites with higher FDR are not used  (0.05)
    -m MIN_WIDTH        The minimum width of hotspots       (50)
__EOF__
  exit 2
}

# Input parsing
MIN_HOTSPOT_WIDTH=50
HOTSPOT_FDR_THRESHOLD=0.05
CUTCOUNTS=

AWK_EXE=$(which mawk 2>/dev/null || which awk)

while getopts 'hc:f:m:' opt ; do
  case "$opt" in
    h) usage ;;
    c) CUTCOUNTS=$OPTARG ;;
    f) HOTSPOT_FDR_THRESHOLD=$OPTARG ;;
    m) MIN_HOTSPOT_WIDTH=$OPTARG ;;
  esac
done

shift $((OPTIND - 1 ))

if [[ $# -lt 2 || -z "$CUTCOUNTS" ]] ; then
  usage
fi

infile=$1
outfile=$2


# Function definitions

# We choose to cap all P-values at 1e-100, i.e. -log10(P) = 100.
# To combat issues awk sometimes has with numbers below 1e-300,
# we parse the FDR initially as a string, and assume anything below 1e-100
# passes the user's threshold.
filter(){
  "$AWK_EXE" \
    -v "threshold=$HOTSPOT_FDR_THRESHOLD" \
    '{
       split($6,y,"-")
       if ((length(y)!=1 && y[2]>100) || $6 <= threshold ){
         print $1"\t"$2"\t"$3
       }
     }'
}

merge(){
  bedops -m - \
    | "$AWK_EXE" \
      -v "minW=$MIN_HOTSPOT_WIDTH" \
      'BEGIN{
         FS="\t"
         OFS="\t"
       }
       {
         chrR=$1
         begPosR=$2
         endPosR=$3
         widthR=endPosR-begPosR

         if(NR>1) {
           if (chrR == chrL) {
             distLR = begPosR - endPosL;
           } else {
             distLR=999999999;
           }
           if (distLR > minW) {
             if (widthL >= minW) {
               print chrL, begPosL, endPosL;
             }
           } else {
             if (widthL < minW || widthR < minW) {
               begPosR = begPosL;
               widthR = endPosR - begPosR;
             } else {
               print chrL, begPosL, endPosL;
             }
           }
         }
         chrL = chrR;
         begPosL = begPosR;
         endPosL = endPosR;
         widthL = widthR;
       }
       END{
         if(widthL >= minW){
           print chrL, begPosL, endPosL
         }
       }'
}

# We report the largest -log10(P) observed at any bp of a hotspot
# as the "score" of that hotspot, where P is the site-specific P-value.
# P-values of 0 will be encountered, and we don't want to do log(0).
# Nonzero P-values as low as 1e-308 have been seen during testing.
# We choose to cap all P-values at 1e-100, i.e. -log10(P) = 100.
# The constant c below converts from natural logarithm to log10.
annotate(){
  "$AWK_EXE" \
    'BEGIN{
      OFS="\t"
      c=-0.4342944819
      max=100
    }
    {
      if (0 == $4) {
        print $1, $2, $3, "id-"NR, $5, ".","-1","-1", max
      } else {
        split($4,y,"-");
        if(length(y)!=1 && y[2] > max) {
          print $1, $2, $3, "id-"NR, $5, ".","-1","-1", max
        } else {
          print $1, $2, $3, "id-"NR, $5, ".","-1","-1", c*log($4)
        }
      }
    }'
}

# Main

unstarch "$infile" \
    | filter \
    | merge \
    | bedmap --faster --sweep-all --delim "\t" --echo --min   - "$infile" \
    | bedmap --faster --sweep-all --delim "\t" --echo --count - "$CUTCOUNTS" \
    | annotate \
    | starch - \
> "$outfile"
