#!/bin/bash -e

if [ $# -lt 2 ] || [ $# -gt 5 ]; then
    echo "usage v1: dnase_index_bwa.sh <genome> <reference.fasta.gz> [<skip_indexing> <mappable_only.starch> [<blacklist.bed.gz>]]"
    echo "Indexes reference for bwa alignment.  Optionally creates hotspot mappable regions. Is independent of DX and encodeD."
    echo "Requires bwa on path.  Making mappable regions needs: faSize, extractCenterSites.sh,"
    echo "bedops (bedmap,sort-bed,starch,starchcat,unstarch) on path."
    exit -1; 
fi
genome=$1               # Genome assembly (e.g. "GRCh38").
ref_fa_gz=$2            # Reference fasta file (e.g. GRCh38.fa.gz).  Will be uncompressed if not already.
skip_indexing="false"
mappable_only_starch="none"
blacklist_bed_gz="none"
if [ $# -gt 3 ]; then
    skip_indexing=$3            # If making mappable targets, optionally skip indexing
    mappable_only_starch=$4     # Mappable regions only file (e.g. GRCh38_no_alts.K36.mappable_only.starch).
    if [ $# -gt 4 ]; then
        blacklist_bed_gz=$5     # Assembly blacklist file (e.g. wgEncodeDacMapabilityConsensusExcludable.bed.gz).
    fi
fi

index_file="${genome}_bwa_index.tgz"
echo "-- Index file will be: '$index_file'"

ref_fa=${ref_fa_gz%.gz}
if [ "$ref_fa" != "$ref_fa_gz" ]; then
    echo "-- Uncompressing reference"
    set -x
    gunzip $ref_fa_gz 
    set +x
fi
if [ "$ref_fa" != "$genome.fa" ]; then
    set -x
    mv $ref_fa ${genome}.fa
    set +x
fi

# Optionally create a mappable regions tar: faSize, sort-bed bedops starch unstarch extractCenterSites.sh
mappable_tar=""
if [ -f $mappable_only_starch ]; then
    mappable_tar="${genome}_hotspot2_v2.0_mappable.tgz"
    echo "-- Create hotspot2 mappable regions archive: ${mappable_tar}"
    echo "-- Generating chrom_sizes.bed from fasta file..."
    set -x
    faSize -detailed ${genome}.fa | awk '{printf "%s\t0\t%s\n",$1,$2}' | sort-bed - > chrom_sizes.bed
    set +x
    #cat $chrom_sizes | awk '{printf "%s\t0\t%s\n",$1,$2}' | sort-bed - > chrom_sizes.bed

    blacklist_bed=""
    if [ -f $blacklist_bed_gz ]; then
        blacklist_bed=${blacklist_bed_gz%.gz}
        if [ "$blacklist_bed" != "$blacklist_bed_gz" ]; then
            echo "-- Uncompressing blacklist..."
            set -x
            gunzip $blacklist_bed_gz 
            set +x
        fi
        echo "-- Subtracting blacklist from mappable regions..."
        set -x
        bedops --difference $mappable_only_starch $blacklist_bed | starch - > mappable_target.starch
        set +x
    else
        echo "-- Using supplied mappable regions as mappable target file..."
        set -x
        cp $mappable_only_starch mappable_target.starch
        set +x
    fi

    echo "-- Creating centerSites file..."
    set -x
    extractCenterSites.sh -c chrom_sizes.bed -M mappable_target.starch -o center_sites.starch
    set +x

    echo "-- Archiving results"
    set -x
    tar -czf $mappable_tar center_sites.starch mappable_target.starch chrom_sizes.bed $mappable_only_starch $blacklist_bed
    set +x
fi

if [ "$skip_indexing" != "true" ]; then
    echo "-- Build index..."
    set -x
    bwa index -p $genome -a bwtsw ${genome}.fa
    set +x

    echo "-- tar and gzip index..."
    set -x
    tar -czf $index_file ${genome}.*
    set +x
        
    if [ "$ref_fa" != "${genome}.fa" ]; then
        set -x
        mv ${genome}.fa $ref_fa
        set +x
    fi
fi

echo "-- The results..."
if [ -f $index_file ]; then
    ls -l $index_file
fi
if [ -f $mappable_tar ]; then
    ls -l $mappable_tar
fi
df -k .

