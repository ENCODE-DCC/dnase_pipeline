#!/usr/bin/env python2.7
# versions.py 0.0.1
#
# Creates versions json string for a particular applet

import sys, os, argparse, json, commands

# APP_TOOLS is a dict keyed by applet script name with a list of tools that it uses.
APP_TOOLS = { 
    "index-bwa.sh":     [ "bwa", "samtools" ],
    "align-bwa-pe.sh":  [ "bwa", "samtools" ],
    "align-bwa-se.sh":  [ "bwa", "samtools" ],
    "bam-filter-pe.sh": [
                            "samtools","edwBamFilter","edwBamStats",#"R",
                            "Rscript","phantompeakqualtools","caTools","snow","spp","gawk","bedtools"
                        ],
    "bam-filter-se.sh": [
                            "samtools","edwBamFilter","edwBamStats",#"R",
                            "Rscript","phantompeakqualtools","caTools","snow","spp","gawk","bedtools"
                        ],
    "call-hotspots.sh": [
                            "hotspot","hotspot.py (GCAP)","samtools",
                            "bedops","bedmap (bedops)","sort-bed (bedops)",
                            "starch (bedops)","starchcat (bedops)","unstarch (bedops)",
                            "bedtools","bamToBed (bedtools)","intersectBed (bedtools)","shuffleBed (bedtools)",
                            "bedToBigBed","bedGraphToBigWig","bedGraphPack"
                        ],
    "sample-hotspots.sh": [
                            "edwBamStats","hotspot","hotspot.py (GCAP)","samtools",
                            "bedops","bedmap (bedops)","sort-bed (bedops)",
                            "starch (bedops)","starchcat (bedops)","unstarch (bedops)",
                            "bedtools","bamToBed (bedtools)","shuffleBed (bedtools)"
                        ],
    "merge-replicates.sh": [ "samtools","bedtools","bigBedToBed","bedToBigBed","bigWigCorrelate","edwComparePeaks" ],
    "fastq-stats.sh": [ "fastqStatsAndSubsample" ]
 }

# ALL_TOOLS contains the printable tool name (key) and the command that is used to determine the version.
ALL_TOOLS = { 
            "bamToBed (bedtools)":      "bamToBed -h 2>&1 | grep Version | awk '{print $2}'",
            "bedGraphPack":             "bedGraphPack 2>&1 | grep 'bedGraphPack v' | awk '{print $2}'",
            "bedGraphToBigWig":         "bedGraphToBigWig 2>&1 | grep 'bedGraphToBigWig v' | awk '{print $2$3}'",
            "bedmap (bedops)":          "bedmap --version 2>&1 | grep version | awk '{print $2}'",
            "bedops":                   "bedops --version 2>&1 | grep version | awk '{print $2}'",
            "bedToBigBed":              "bedToBigBed 2>&1 | grep 'bedToBigBed v' | awk '{print $2$3}'",
            "bedtools":                 "bedtools --version 2>&1 | awk '{print $2}'",
            "bigBedToBed":              "bigBedToBed 2>&1 | grep 'bigBedToBed v' | awk '{print $2}'",
            "bigWigCorrelate":          "echo unversioned", #"bigWigCorrelate 2>&1 | grep 'bigWigCorrelate -' | awk '{print $3,$4,$5}'",
            "bwa":                      "bwa 2>&1 | grep Version | awk '{print $2}'",
            "caTools":                  "grep caTools_ phantompeakqualtools/install.log | head -1 | sed 's/_/ /' | awk '{print $4}' | sed 's/\.tar\.gz.*//'",
            "edwBamFilter":             "edwBamFilter 2>&1 | grep 'edwBamFilter v' | awk '{print $2}'",
            "edwBamStats":              "edwBamStats 2>&1 | grep 'edwBamStats v' | awk '{print $2}'",
            "edwComparePeaks":          "echo unversioned", #"edwComparePeaks 2>&1 | grep 'edwComparePeaks -' | awk '{print $3,$4,$5,$6}'",
            "fastqStatsAndSubsample":   "fastqStatsAndSubsample 2>&1 | grep 'fastqStatsAndSubsample v' | awk '{print $2}'",
            "gawk":                     "gawk --version | grep Awk | awk '{print $3}'",
            "hotspot":                  "hotspot 2>&1 | grep HotSpot | awk '{print $1}'",
            "hotspot.py (GCAP)":        "python2.7 /usr/bin/hotspot.py -h | grep Version | awk '{print $8}'",
            "intersectBed (bedtools)":  "intersectBed 2>&1 | grep Version | awk '{print $2}'",
            "phantompeakqualtools":     "grep Version phantompeakqualtools/README.txt | awk '{print $2}'",
            "R":                        "R --version | grep 'R version' | awk '{print $3,$4}'",
            "Rscript":                  "Rscript --version 2>&1 | awk '{print $5,$6}'",
            "samtools":                 "samtools 2>&1 | grep Version | awk '{print $2}'",
            "shuffleBed (bedtools)":    "shuffleBed -h 2>&1 | grep Version | awk '{print $2}'",
            "snow":                     "grep snow_ phantompeakqualtools/install.log | head -1 | sed 's/_/ /' | awk '{print $4}' | sed 's/\.tar\.gz.*//'",
            "spp":                      "grep spp_ phantompeakqualtools/installPkgs.R | sed 's/_/ /' | awk '{print $2}' | sed 's/\.tar\.gz.*//'",
            "sort-bed (bedops)":        "sort-bed --version 2>&1 | grep version | awk '{print $2}'",
            "starch (bedops)":          "starch --version 2>&1 | grep version | awk '{print $3}'",
            "starchcat (bedops)":       "starchcat --version 2>&1 | grep version | awk '{print $3}'",
            "unstarch (bedops)":        "unstarch --version 2>&1 | grep version | awk '{print $3}'"
            }

def main():
    parser = argparse.ArgumentParser(description =  "Versions parser for a dx applet. " + \
                                                    "Prints version lines to stderr and json string to stdout.")
    parser.add_argument('-a','--applet', required=True,
                        help="Applet to print versions for")
    parser.add_argument('-av','--appver', required=True,
                        help="Version of applet")
    parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False, 
                        help="Don't print versions to stderr.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False, 
                        help="Show the command-line that is used to get the version.")


    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 3:
        parser.print_usage()
        return
        
    versions = {}
    versions["DX applet"] = { args.applet: args.appver }
    if not args.quiet:
        sys.stderr.write("********\n")
        sys.stderr.write("* Running " + args.applet + ": " + args.appver+ "\n")
        
    tools = APP_TOOLS[args.applet]
    for tool in tools:
        cmd = ALL_TOOLS[tool]
        if args.verbose:
            sys.stderr.write("cmd> " + cmd + "\n")
        err, ver = commands.getstatusoutput(cmd)
        versions[tool] = ver
        if not args.quiet:
            sys.stderr.write("* " + tool + " version: " + ver + "\n")

    if not args.quiet:
        sys.stderr.write("********\n")
    
    print json.dumps(versions) 
     
if __name__ == '__main__':
    main()


