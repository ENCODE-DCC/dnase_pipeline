#!/bin/bash
# sample-hotspots

script_name="sample-hotspots.sh"
script_ver="0.2.0"

main() {
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
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of bam_to_sample: '$bam_to_sample'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of read_length: '$read_length'"
    echo "* Value of genome: '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_to_sample" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_to_sample" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    sample_root="${bam_root}_sample_5M"
    echo "* out: sample file: '${sample_root}.bam'"

    dx download "$chrom_sizes" -o chromSizes.txt
    # sort-bed is important!
    grep -v chrM chromSizes.txt | awk '{printf "%s\t0\t%s\t%s\n",$1,$2,$1}' | sort-bed - > ${genome}.chromInfo.bed

    read_size=`dx describe "$bam_to_sample" --details --json | grep "\"average length\":" | awk '{print $3}' | tr -d ,`
    if [ "$read_size" == "" ]; then
        read_size=`dx describe "$bam_to_sample" --details --json | grep "\"readSizeMean\":" | awk '{print $2}' | tr -d ,`
        if [ "$read_size" == "" ] && [ -f /usr/bin/parse_property.py ]; then
            read_size=`parse_property.py -f "$bam_to_sample" -k "average length" --quiet`
            if [ "$read_size" == "" ]; then
                read_size=`parse_property.py -f "$bam_to_sample" -s edwBamStats -k readSizeMean --quiet`
            fi
        fi
    fi
    if [ "$read_size" != "" ] && [ "$read_length" -ne "$read_size" ]; then
        if [ "$read_size" == "32" ] || [ "$read_size" == "36" ] || [ "$read_size" == "40" ] || [ "$read_size" == "50" ] \
        || [ "$read_size" == "58" ] || [ "$read_size" == "72" ] || [ "$read_size" == "76" ] || [ "$read_size" == "100" ]; then
            echo "* NOTE: Read length ($read_length) does not match discovered read size ($read_size). Using $read_size."
            read_length=$read_size
        else
            echo "* WARNING: Read length ($read_length) does not match discovered read size ($read_size)."
        fi
    fi

    # TODO: Need to make bam.bai?
    echo "* Indexing bam..."
    set -x
    samtools index ${bam_root}.bam
    set +x
    
    echo "* Sampling bam..."
    set -x
    edwBamStats -sampleBamSize=5000000 -u4mSize=5000000 -sampleBam=${sample_root}.bam ${bam_root}.bam ${sample_root}_stats.txt
    samtools index ${sample_root}.bam
    set +x
    
    echo "* Running hotspot.py..."
    set -x
    mkdir tmp
    mkdir out
    cp ${genome}.chromInfo.bed hotspot/hotspot-distr/data/
    # May also need to do something about "Satellite.${genome}.bed"
    #cp /usr/bin/Satellite.${genome}.bed hotspot/hotspot-distr/data   # hg19 version already there!
    mappable=${genome}.K${read_length}.mappable_only
    wget http://www.uwencode.org/proj/hotspot/${mappable}.bed -O hotspot/hotspot-distr/data/${mappable}.bed >> install.log 2>&1
    python2.7 /usr/bin/hotspot.py -o hotspot/hotspot-distr/ ${sample_root}.bam $genome DNase-seq $read_length tmp out
    # hanging on: /home/dnanexus/hotspot/hotspot-distr/pipeline-scripts/run_generate_random_lib

    cp tmp/${sample_root}.spot.out ${sample_root}_hotspot_qc.txt
    set +x

    echo "* Prepare metadata..."
    qc_sampled=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_sampled=`qc_metrics.py -n edwBamStats -f ${sample_root}_stats.txt`
        meta=`qc_metrics.py -n hotspot -f ${sample_root}_hotspot_qc.txt`
        qc_sampled=`echo $qc_sampled, $meta`
    fi

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_sample_5M=$(dx upload ${sample_root}.bam --details "{ $qc_sampled }" --property QC="{ $qc_sampled }" --property SW="$versions" --brief)
    hotspot_sample_5M_qc=$(dx upload ${sample_root}_hotspot_qc.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_sample_5M "$bam_sample_5M" --class=file
    dx-jobutil-add-output hotspot_sample_5M_qc "$hotspot_sample_5M_qc" --class=file
    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
