#!/bin/env python3

import sys
import gzip

def transform_line(line):

    line = line.strip()

    if line.startswith("@"):

        insert_loc = line.find(' ')
        umi_loc = line.find('+') + 1
        umi = line[umi_loc:] if umi_loc else ""

        line = "%s#%s%s" % (line[:insert_loc], umi, line[insert_loc:])

    return line

def transform_file(infilename, outfilename):

    outfile = gzip.open(outfilename, 'wt')
    infile = gzip.open(infilename, 'rt')

    for line in infile:
        outfile.write( transform_line(line) )
        outfile.write('\n')

    outfile.close()
    infile.close()

if __name__ == "__main__":

    transform_file(sys.argv[1], sys.argv[2])
