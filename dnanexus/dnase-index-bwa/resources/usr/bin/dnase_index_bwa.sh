#!/bin/bash -e

if [ $# -lt 2 ] ||  [ $# -gt 3 ]; then
    echo "usage v1: dnase_index_bwa.sh <reference.fasta.gz> <genome> [<gender>]"
    echo "Indexes reference for bwa alignment.  Is independent of DX and encodeD."
    echo "Requires bwa on path."
    exit -1; 
fi
ref_fa_gz=$1  # Reference fasta file (e.g. GRCh38.fa.gz).  Will be uncompressed if not already.
genome=$2     # Genome assembly (e.g. "GRCh38").
if [ $# -eq 3 ]; then
    gender=$3     # Gender of assembly (e.g. "female").  Anything besides "female" or "XX" is interpreted as "XY".
    if [ "$gender" != "female" ] && [ "$gender" != "XX" ]; then
        gender="XY"
    fi
    index_id="${genome}_${gender}"
else
    index_id=$genome
fi

index_file="${index_id}_bwa_index.tgz"
echo "-- Index file will be: '$index_file'"

ref_fa=${ref_fa_gz%.gz}
if [ "$ref_fa" != "$ref_fa_gz" ]; then
    echo "-- Uncompressing reference"
    set -x
    gunzip $ref_fa_gz 
    set +x
fi
if [ "$ref_fa" != "$index_id.fa" ]; then
    set -x
    mv $ref_fa ${index_id}.fa
    set +x
fi

echo "-- Build index..."
set -x
bwa index -p $index_id -a bwtsw ${index_id}.fa
set +x

echo "-- tar and gzip index..."
set -x
tar -czf $index_file ${index_id}.*
set +x
    
if [ "$ref_fa" != "$index_id.fa" ]; then
    set -x
    mv ${index_id}.fa $ref_fa
    set +x
fi

#if [ "$ref_fa" != "$ref_fa_gz" ]; then
#    echo "-- Recompressing reference"
#    set -x
#    gzip $ref_fa 
#    set +x
#fi

echo "-- The results..."
ls -l $index_file
df -k .

