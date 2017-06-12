#!/usr/bin/env python2.7
# tool_versions.py v1.1  Creates "SW" versions json string for a particular DX applet.
#                        Write request to stdout and verbose info to stderr.  This allows easy use
#                        in dx app scripts.
#
# Creates versions json string for a particular applet

import sys
import argparse
import json
import commands

# APP_TOOLS is a dict keyed by applet script name with a list of tools that it uses.
APP_TOOLS = {
    "dnase-index-bwa":    ["dnase_index_bwa.sh", "bwa", "hotspot2", "bedops"],
    "dnase-align-bwa-pe": ["dnase_align_bwa_pe.sh", "bwa", "samtools", "edwBamStats",
                           "trim-adapters-illumina", "fastq_umi_add.py (stampipes)"],
    "dnase-align-bwa-se": ["dnase_align_bwa_se.sh", "bwa", "samtools", "edwBamStats",
                           "cutadapt"],
    "dnase-filter-pe":    ["dnase_filter_pe.sh", "samtools", "filter_reads.py (stampipes)",
                           "java", "picard"],
    "dnase-filter-se":    ["dnase_filter_se.sh", "samtools", "java", "picard"],
    "dnase-qc-bam":       ["dnase_qc_bam.sh", "samtools", "edwBamFilter", "edwBamStats",  # "R",
                           "Rscript", "phantompeakqualtools", "caTools", "snow", "spp", "gawk",
                           "hotspot1", "hotspot.py", "bedops", "bedtools"],
    "dnase-density":      ["dnase_density.sh", "samtools", "bedops", "bedGraphToBigWig", "gawk"],
    "dnase-call-hotspots": ["dnase_hotspot.sh", "samtools", "hotspot2", "bedops", "modwt", "gawk", "mawk",
                            "bedToBigBed", "bedGraphToBigWig"],
    "dnase-rep-corr":      ["dnase_rep_corr.sh", "chromCor.Rscript", "bigWigToWig", "bedops"],
    # special for optional output
    "dnase-index-bwa(hotspot2)": ["dnase_index_bwa.sh", "hotspot2", "bedops"],
    }

# Virtual apps only differ from their parent by name/version. 
VIRTUAL_APPS = {
    "dnase-qc-bam-alt":          "dnase-qc-bam",
    "dnase-density-alt":         "dnase-density",
    "dnase-call-hotspots-alt":   "dnase-call-hotspots",
    "dnase-idr-alt":             "dnase-idr",
    "dnase-rep-corr-alt":        "dnase-rep-corr",
    }

# ALL_TOOLS contains printable tool name (key) and the command that is used to determine version.
ALL_TOOLS = {"Anaconda3":                "ls Anaconda3*.sh | head -1 | cut -d - -f 2",
             "bedGraphPack":             "bedGraphPack 2>&1 | grep 'bedGraphPack v' | awk '{print $2}'",
             "bedGraphToBigWig":         "bedGraphToBigWig 2>&1 | grep 'bedGraphToBigWig v' | awk '{print $2$3}'",
             "bedops":                   "bedops --version 2>&1 | grep version | awk '{print $2}'",
             # "bam2bed (bedops)":         "bedops --version 2>&1 | grep version | awk '{print $2}'", # Note: no version.. subsituting bedops
             # "bedmap (bedops)":          "bedmap --version 2>&1 | grep version | awk '{print $2}'",
             # "convert2bed (bedops)":     "convert2bed --version 2>&1 | grep version | awk '{print $2}'",
             # "sort-bed (bedops)":        "sort-bed --version 2>&1 | grep version | awk '{print $2}'",
             # "starch (bedops)":          "starch --version 2>&1 | grep version | awk '{print $3}'",
             # "starchcat (bedops)":       "starchcat --version 2>&1 | grep version | awk '{print $3}'",
             # "unstarch (bedops)":        "unstarch --version 2>&1 | grep version | awk '{print $3}'",
             "bedToBigBed":              "bedToBigBed 2>&1 | grep 'bedToBigBed v' | awk '{print $3}'",
             "bedtools":                 "bedtools --version 2>&1 | awk '{print $2}'",
             # "bamToBed (bedtools)":      "bamToBed -h 2>&1 | grep Version | awk '{print $2}'",
             # "intersectBed (bedtools)":  "intersectBed 2>&1 | grep Version | awk '{print $2}'",
             # "shuffleBed (bedtools)":    "shuffleBed -h 2>&1 | grep Version | awk '{print $2}'",
             "bigBedToBed":              "bigBedToBed 2>&1 | grep 'bigBedToBed v' | awk '{print $2}'",
             "bigWigCorrelate":          "md5sum /usr/bin/bigWigCorrelate | awk '{printf \"unversioned %-8.8s\",$1}'",
             "bwa":                      "bwa 2>&1 | grep Version | awk '{print $2}'",
             "caTools":                  "grep caTools_ phantompeakqualtools/install.log | head -1 | sed 's/_/ /' | awk '{print $4}' | sed 's/\.tar\.gz.*//'",
             "edwBamFilter":             "edwBamFilter 2>&1 | grep 'edwBamFilter v' | awk '{print $2}'",
             "edwBamStats":              "edwBamStats 2>&1 | grep 'edwBamStats v' | awk '{print $2}'",
             "edwComparePeaks":          "md5sum /usr/bin/edwComparePeaks | awk '{printf \"unversioned %-8.8s\",$1}'", #"edwComparePeaks 2>&1 | grep 'edwComparePeaks -' | awk '{print $3,$4,$5,$6}'",
             "faSize":                   "md5sum /usr/bin/faSize | awk '{printf \"unversioned %-8.8s\",$1}'", #"bigWigCorrelate 2>&1 | grep 'bigWigCorrelate -' | awk '{print $3,$4,$5}'",
             "fastqStatsAndSubsample":   "fastqStatsAndSubsample 2>&1 | grep 'fastqStatsAndSubsample v' | awk '{print $2}'",
             # stampipes:
             "fastq_umi_add.py (stampipes)":         "md5sum /usr/bin/fastq_umi_add.py | awk '{printf \"unversioned %-8.8s\",$1}'", # From https://github.com/StamLab/stampipes/tree/encode-release/scripts/umi/
             "filter_reads.py (stampipes)":          "md5sum ./filter_reads.py | awk '{printf \"unversioned %-8.8s\",$1}'",  # From https://github.com/StamLab/stampipes/tree/encode-release/scripts/bwa/
             # "umi_sort_sam_annotate.awk (stampipes)":"md5sum /usr/bin/umi_sort_sam_annotate.awk | awk '{printf \"unversioned %-8.8s\",$1}'",  # From https://github.com/StamLab/stampipes/tree/encode-release/scripts/umi/
             # "mark_umi_dups.mk (stampipes)":       "md5sum /usr/bin/mark_duplicates.mk | awk '{printf \"unversioned %-8.8s\",$1}'",  # From https://github.com/StamLab/stampipes/tree/encode-release/makefiles/umi/mark_duplicates.mk (dcc minimal change)
             "gawk":                     "gawk --version | grep Awk | awk '{print $3}'",
             "idr":                      "idr/bin/idr --version 2>&1 | grep IDR | awk '{print $2}'",
             "mawk":                     "mawk -W version 2>&1 | grep mawk | awk '{print $2}'",
             "hotspot1":                 "hotspot 2>&1 | grep HotSpot | awk '{printf \"%s-%s\",$1,$2}'",
             "hotspot.py":               "hotspot.py -h | grep Version | awk '{print $8}'",
             "java":                     "java -version 2>&1 | head -1 | awk '{print $3}' | tr -d '\"'",
             "phantompeakqualtools":     "grep Version phantompeakqualtools/README.txt | awk '{print $2}'",
             "R":                        "R --version | grep 'R version' | awk '{print $3,$4}'",
             "Rscript":                  "Rscript --version 2>&1 | awk '{print $5,$6}'",
             "samtools":                 "samtools 2>&1 | grep Version | awk '{print $2}'",
             "snow":                     "grep snow_ phantompeakqualtools/install.log | head -1 | sed 's/_/ /' | awk '{print $4}' | sed 's/\.tar\.gz.*//'",
             "spp":                      "grep spp_ phantompeakqualtools/installPkgs.R | sed 's/_/ /' | awk '{print $2}' | sed 's/\.tar\.gz.*//'",
             #"hotspot2":                 "hotspot2 --version | awk '{print $3}'",
             #"hotspot2":                 "[ -e /usr/bin/hotspot2.version ] && cat /usr/bin/hotspot2.version || hotspot2 --version | awk '{print $3}'",
             "hotspot2":                 "hotspot2_part1 --version | awk '{print $3}'",
             "modwt":                    "md5sum /usr/bin/modwt | awk '{printf \"unversioned %-8.8s\",$1}'", # From https://github.com/StamLab/modwt/tree/1.0
             "picard":                   "java -jar ./picard.jar MarkDuplicates --version", # From https://github.com/broadinstitute/picard.git
             "pigz":                     "pigz --version 2>&1 | awk '{print $2}'",
             "trim-adapters-illumina":   "trim-adapters-illumina --version 2>&1 | awk '{print $3}'", # https://bitbucket.org/jvierstra/bio-tools/get/master.tar.gz https://bitbucket.org/jvierstra/bio-tools/src/6fe54fa5a3d9b5c930ee77e8ccd757b347c86ac1/apps/trim-adapters-illumina/?at=master
             "chromCor.Rscript":         "md5sum /usr/bin/chromCor.Rscript | awk '{printf \"unversioned %-8.8s\",$1}'", # emailed from Richard Sandstrom  Will reside in our github
             "bigWigToWig":              "md5sum /usr/bin/bigWigToWig | awk '{printf \"unversioned %-8.8s\",$1}'",
             "dnase_index_bwa.sh":       "dnase_index_bwa.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_align_bwa_pe.sh":    "dnase_align_bwa_pe.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_align_bwa_se.sh":    "dnase_align_bwa_se.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_filter_pe.sh":       "dnase_filter_pe.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_filter_se.sh":       "dnase_filter_se.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_qc_bam.sh":          "dnase_qc_bam.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_density.sh":         "dnase_density.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_hotspot.sh":         "dnase_hotspot.sh | grep usage | awk '{print $2}' | tr -d :",
             "dnase_rep_corr.sh":        "dnase_rep_corr.sh | grep usage | awk '{print $2}' | tr -d :",
             "cutadapt":                     "cutadapt --version",
            }

def parse_dxjson(dxjson):
    '''Parses the dnanexus-executable.json file in the job directory to get applet name and version.'''
    with open(dxjson) as data_file:
        dxapp = json.load(data_file)

    appver = "unknown"
    applet = dxapp.get("name")
    if "version" in dxapp:
        appver = dxapp.get("version")
    else:
        title = dxapp.get("title")
        last_word = title.split(' ')[-1]
        if last_word.startswith('(virtual-') and last_word.endswith(')'):
            appver = last_word[9:-1]
        elif last_word.startswith('(v') and last_word.endswith(')'):
            appver = last_word[2:-1]

    return (applet, appver)


def main():
    parser = argparse.ArgumentParser(description="Versions parser for a dx applet. " + \
                                                 "Prints version lines to stderr and json string to stdout. " + \
                                                 "MUST specify either --applet and --appver or --dxjson.")
    parser.add_argument('-a', '--applet', required=False,
                        help="Applet to print versions for")
    parser.add_argument('-av', '--appver', required=False,
                        help="Version of applet")
    parser.add_argument('-j', '--dxjson', required=False,
                        help="Use dnanexus json file to discover 'applet' and 'appver'")
    parser.add_argument('-k', '--key',
                        help='Prints just the value for this key.',
                        default=None,
                        required=False)
    parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False,
                        help="Don't print versions to stderr.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False,
                        help="Show the command-line that is used to get the version.")

    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 3:
        parser.print_usage()
        return

    if (args.applet == None or args.appver == None) and args.dxjson == None:
        parser.print_help()
        return

    applet = args.applet
    appver = args.appver

    if args.dxjson != None:
        (applet, appver) = parse_dxjson(args.dxjson)

    versions = {}
    versions["DX applet"] = {applet: appver}
    if not args.quiet:
        sys.stderr.write("********\n")
        sys.stderr.write("* Running " + applet + ": " + appver + "\n")

    if applet in VIRTUAL_APPS:
        tools = APP_TOOLS[VIRTUAL_APPS[applet]]
    else:
        tools = APP_TOOLS[applet]
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

    if args.key != None:
        if args.key in versions:
            print versions[args.key]
            if not args.quiet:
                sys.stderr.write(versions[args.key] + '\n')
        elif args.key in versions["DX applet"]:
            print versions["DX applet"][args.key]
            if not args.quiet:
                sys.stderr.write(versions["DX applet"][args.key] + '\n')
        else:
            print ''
            if not args.quiet:
                sys.stderr.write('(not found)\n')
    else:
        print json.dumps(versions)


if __name__ == '__main__':
    main()


