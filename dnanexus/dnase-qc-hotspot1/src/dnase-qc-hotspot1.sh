#!/bin/bash
# dnase-qc-hotspot1.sh - Calls hotspot1 on a sample for qc for the ENCODE DNase-seq pipeline.

main() {
    # Questions for reviving hotspot1 for SPOT score:
    # Do we use the _sample.bam created by dnase_eval_bam? Yes.
    # Do we try merge bams from 2 replicates? No.
    # How long will it take?  15M sample: 1h16m $1.50
    echo "* Installing hotspot and dependencies (gsl)..." 2>&1 | tee -a install.log
    exe_dir="`pwd`"
    set -x
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz -O gsl.tgz >> install.log 2>&1
    mkdir gsl
    tar -xzf gsl.tgz -C gsl --strip-components=1
    cd gsl
    ./configure >> install.log 2>&1
    make > install.log 2>&1
    sudo make install >> install.log 2>&1
    gsl-config --libs > install.log 2>&1
    export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
    cd ..
    wget https://github.com/rthurman/hotspot/archive/v4.1.0.tar.gz -O hotspot.tgz >> install.log 2>&1
    mkdir hotspot
    tar -xzf hotspot.tgz -C hotspot --strip-components=1
    cd hotspot/hotspot-distr/hotspot-deploy
    make >> install.log 2>&1
    # Can either put bin in path or copy contents of bin to /usr/bin
    export PATH=${exe_dir}/hotspot/hotspot-distr/hotspot-deploy/bin:${PATH}
    cd ../../../
    set +x; 
    # additional executables in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_file:         '$bam_file'"
    echo "* Value of hotspot_mappable: '$hotspot_mappable'"
    #echo "* Value of read_length:      '$read_length'"
    echo "* Value of genome:           '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_file" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_file" -o ${bam_root}.bam
    
    #dx download "$chrom_sizes" -o chrom.sizes
    mappable_archive=`dx describe "$hotspot_mappable" --name`
    dx download "$hotspot_mappable"

    # check read_length
    read_length=36   # Only supporting read_length 36 at this time
    #if [ -f /usr/bin/parse_property.py ]; then 
    #    read_len=`parse_property.py -f "$bam_file" -p "read_length" --quiet`
    #    if [ "$read_len" != "" ] && [ "$read_length" -ne "$read_len" ]; then
    #        if [ "$read_len" == "32" ] || [ "$read_len" == "36" ] || [ "$read_len" == "40" ] || [ "$read_len" == "50" ] \
    #        || [ "$read_len" == "58" ] || [ "$read_len" == "72" ] || [ "$read_len" == "76" ] || [ "$read_len" == "100" ]; then
    #            echo "* NOTE: Read length ($read_length) does not match discovered read size ($read_len). Using $read_len."
    #            read_length=$read_len
    #        else
    #            echo "* WARNING: Read length ($read_length) does not match discovered read size ($read_len)."
    #        fi
    #    fi
    #fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_qc_hotspot1.sh ${bam_root}.bam $genome $mappable_archive $read_length hotspot/
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Prepare metadata..."
    qc_hotspot1=''
    reads_sample=5000000
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_hotspot1=`qc_metrics.py -n hotspot -f ${bam_root}_hotspot1_qc.txt`
    fi

    echo "* Upload results..."
    bam_hotspot1_qc=$(dx upload ${bam_root}_hotspot1_qc.txt --details "{ $qc_hotspot1 }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_hotspot1_qc "$bam_hotspot1_qc" --class=file

    dx-jobutil-add-output metadata "{ $qc_hotspot1 }" --class=string

    echo "* Finished."
}
