#!/bin/bash
# index-bwa.sh

script_name="index-bwa.sh"
script_ver="0.1.0"

main() {
    # Executable in resources/usr/bin

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi
    #echo "*****"
    #echo "* Running: index-bwa.sh [v0.1.0]"
    #echo "* bwa version: "`bwa 2>&1 | grep Version | awk '{print $2}'`
    #echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    #echo "*****"

    echo "* Value of reference: '$reference'"
    echo "* Value of genome: '$genome'"
    echo "* Value of gender: '$gender'"

    index_id="${genome}_${gender}"
    index_file="${index_id}_bwa_index.tgz"
    echo "* Index file will be: '$index_file'"

    echo "* Download files..."
    #ref_root=`dx describe "$reference" --name`
    #ref_root=${ref_root%.fasta.gz}
    #ref_root=${ref_root%.fa.gz}
    dx download "$reference" -o "$index_id".fa.gz
    gunzip "$index_id".fa.gz 
    ref="$index_id".fa

    echo "* Reference file: '$ref'"

    # Fill in your application code here.

    echo "* Build index..."
    set -x
    bwa index -p $index_id -a bwtsw $ref
    ls 
    ls -l ${index_id}.*
    set +x
    
    echo "* tar and gzip index..."
    set -x
    tar -czf $index_file ${index_id}.*
    set +x
    
    echo "* Upload Results..."
    bwa_version=`bwa 2>&1 | grep Version | awk '{print $2}'`
    bwa_index=$(dx upload $index_file --property genome="$genome" --property gender="$gender" --property SW="$versions" --brief)

    dx-jobutil-add-output bwa_index $bwa_index --class=file
    echo "* Finished."
}
