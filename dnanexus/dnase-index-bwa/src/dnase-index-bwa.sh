#!/bin/bash
# dnase-index-bwa.sh

main() {
    # Executable in resources/usr/bin

    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
        appver=`tool_versions.py -q --dxjson dnanexus-executable.json --key "dnase-index-bwa"`
        mappable_versions=`tool_versions.py --applet "dnase-index-bwa(hotspot2)" --appver $appver`
    fi

    echo "* Value of reference:        '$reference'"
    echo "* Value of mappable_regions: '$mappable_regions'"
    echo "* Value of blacklist:        '$blacklist'"
    echo "* Skip indexing:             '$skip_indexing'"
    # Prefer to discover genome and gender
    source_msg="Value of"
    gender_msg="Default "
    gender="XY" # default 
    if [ -f /usr/bin/parse_property.py ]; then 
        genome_prop=`parse_property.py -f "$reference" -p "genome" --quiet`
        gender_prop=`parse_property.py -f "$reference" -p "gender" --quiet`
        if [ "$gender_prop" != "" ]; then
            genome=$genome_prop
            source_msg="Discovered"
        fi
        if [ "$gender_prop" != "" ]; then
            gender=$gender_prop
            gender_msg="Discovered"
        fi
    fi
    if [ "$genome" == "" ]; then
        echo "Reference genome could not be determined and must be supplied as arguments."
        exit 1
    fi
    echo "* ${source_msg} genome: '$genome'"
    echo "* ${gender_msg} gender: '$gender'"

    echo "* Download files..."
    dx download "$reference" -o ${genome}.fa.gz
    echo "* Reference file: '${genome}.fa.gz'"

    mappable_starch=""
    blacklist_bed_gz=""
    if [ "$mappable_regions" != "" ]; then
        mappable_starch=`dx describe "$mappable_regions" --name`
        dx download "$mappable_regions"
        echo "* Mappable regions file: '$mappable_starch'"
        if [ "$blacklist" != "" ]; then
            blacklist_bed_gz=`dx describe "$blacklist" --name`
            dx download "$blacklist"
            echo "* Blacklist file: '$$blacklist'"
        fi
        mappable_root=${mappable_starch%.bed.gz}
        if [ "$mappable_root" != "$mappable_starch" ]; then
            echo "* Converting $mappable_root.bed.gz to $mappable_root.starch..."
            gunzip $mappable_root.bed.gz
            starch $mappable_root.bed > $mappable_root.starch
            mappable_starch=$mappable_root.starch
        fi
    fi
    
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_index_bwa.sh $genome ${genome}.fa.gz $skip_indexing $mappable_starch $blacklist_bed_gz
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    index_file="${genome}_bwa_index.tgz"
    mappable_file="${genome}_hotspot_mappable.tgz"

    echo "* Upload Results..."
    if [ -f $index_file ]; then
        bwa_index=$(dx upload $index_file --property genome="$genome" --property gender="$gender" --property SW="$versions" --brief)
        dx-jobutil-add-output bwa_index $bwa_index --class=file
    fi

    if [ -f $mappable_file ]; then
        mappable_tar=$(dx upload $mappable_file --property genome="$genome" --property gender="$gender" --property SW="$mappable_versions" --brief)
        dx-jobutil-add-output mappable_tar $mappable_tar --class=file
    fi

    echo "* Finished."
}
