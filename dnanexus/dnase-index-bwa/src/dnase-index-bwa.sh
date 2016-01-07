#!/bin/bash
# dnase-index-bwa.sh

main() {
    # Executable in resources/usr/bin

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of reference: '$reference'"
    # Prefer to discover genome and gender
    source_msg="Value of"
    if [ -f /usr/bin/parse_property.py ]; then 
        genome_prop=`parse_property.py -f "$reference" -p "genome" --quiet`
        gender_prop=`parse_property.py -f "$reference" -p "gender" --quiet`
        if [ "$genome_prop" != "" ] &&  [ "$gender_prop" != "" ]; then
            genome=$genome_prop
            gender=$gender_prop
            source_msg="Discovered"
        fi
    fi
    if [ "$genome" == "" ] || [ "$gender" == "" ]; then
        echo "Reference genome and/or gender could not be determined and must be supplied as arguments."
        exit 1
    fi
    echo "* ${source_msg} genome: '$genome'"
    echo "* ${source_msg} gender: '$gender'"

    index_id="${genome}_${gender}"

    echo "* Download files..."
    dx download "$reference" -o "$index_id".fa.gz

    echo "* Reference file: '$ref'"

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_index_bwa.sh "$index_id".fa.gz $genome $gender
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    index_file="${index_id}_bwa_index.tgz"

    echo "* Upload Results..."
    bwa_index=$(dx upload $index_file --property genome="$genome" --property gender="$gender" --property SW="$versions" --brief)

    dx-jobutil-add-output bwa_index $bwa_index --class=file
    echo "* Finished."
}
