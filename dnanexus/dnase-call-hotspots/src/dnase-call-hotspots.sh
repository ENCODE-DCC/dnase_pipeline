#!/bin/bash
# dnase-call-hotspots.sh - Call peaks with 'hotspot' for the ENCODE DNase-seq pipeline.

script_name="dnase-call-hotspots.sh"
script_ver="0.2.1"

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

    echo "* Value of bam_to_call: '$bam_to_call'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of read_length: '$read_length'"
    echo "* Value of genome:      '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_to_call" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_to_call" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    narrowPeak_root="${bam_root}_narrowPeak_hotspot"
    broadPeak_root="${bam_root}_broadPeak_hotspot"
    signal_root="${bam_root}_signal_hotspot"
    echo "* out: narrowPeak files: '${narrowPeak_root}.bed' and '${narrowPeak_root}.bb'"
    echo "* out: broadPeak files: '${broadPeak_root}.bed' and '${broadPeak_root}.bb'"
    echo "* out: signal file: '${signal_root}.bw'"

    dx download "$chrom_sizes" -o chromSizes.txt
    # sort-bed is important!
    grep -v chrM chromSizes.txt | awk '{printf "%s\t0\t%s\t%s\n",$1,$2,$1}' | sort-bed - > ${genome}.chromInfo.bed

    read_len=`parse_property.py -f "$bam_to_call" -p "read_length" --quiet`
    if [ "$read_len" == "" ]; then
        echo "* Running edwBamStats on '${bam_root}.bam'"
        set -x
        edwBamStats ${bam_root}.bam ${bam_root}_edwBamStats.txt
        set +x
        read_len=`qc_metrics.py -n edwBamStats -f ${bam_root}_edwBamStats.txt -k readSizeMean`
    fi
    if [ "$read_len" != "" ] && [ "$read_length" -ne "$read_len" ]; then
        if [ "$read_len" == "32" ] || [ "$read_len" == "36" ] || [ "$read_len" == "40" ] || [ "$read_len" == "50" ] \
        || [ "$read_len" == "58" ] || [ "$read_len" == "72" ] || [ "$read_len" == "76" ] || [ "$read_len" == "100" ]; then
            echo "* NOTE: Read length ($read_length) does not match discovered read size ($read_len). Using $read_len."
            read_length=$read_len
        else
            echo "* WARNING: Read length ($read_length) does not match discovered read size ($read_len)."
        fi
    fi

    # TODO: Need to make bam.bai?
    echo "* Indexing bam..."
    set -x
    samtools index ${bam_root}.bam
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
    python2.7 /usr/bin/hotspot.py hotspot/hotspot-distr/ ${bam_root}.bam $genome DNase-seq $read_length tmp out
    cp tmp/${bam_root}.spot.out ${bam_root}_hotspot_out.txt
    set +x

    echo "* Generating narrowPeaks..."
    set -x
    # Convert output into ENCODE formats,
    # Convert narrowPeaks.bed and several stray columns to narrowPeaks.bigBed in $4
    paste out/narrowPeaks.bed out/narrowPeaks.dens out/narrowPeaks.pval | \
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "p" NR, 0, ".", $4, $5, -1, -1}' - | \
        sort -k1,1 -k2,2n - > ${narrowPeak_root}.bed
    bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${narrowPeak_root}.bed chromSizes.txt ${narrowPeak_root}.bb
    wc -l ${narrowPeak_root}.bed > ${narrowPeak_root}_qc.txt
    set +x

    echo "* Generating broadPeaks..."
    set -x
    # Convert broadPeaks.bed and broadPeaks.pVal to broadPeaks.bigBed in $5
    paste out/broadPeaks.bed out/broadPeaks.pval | \
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "hot" NR, 0, ".", $5, $6, -1}' - | \
        sort -k1,1 -k2,2n - > ${broadPeak_root}.bed
    bedToBigBed -as=/usr/bin/broadPeak.as -type=bed6+3 ${broadPeak_root}.bed chromSizes.txt ${broadPeak_root}.bb
    wc -l ${broadPeak_root}.bed > ${broadPeak_root}_qc.txt
    set +x

    echo "* Generating signal..."
    set -x
    # Convert starched bedGraph to mappable-only bigWig in $6
    unstarch out/density.bed.starch > tmp/tmp.bed
    intersectBed -a tmp/tmp.bed -b hotspot/hotspot-distr/data/${genome}.K${read_length}.mappable_only.bed -f 1.00 | \
        cut -f 1,2,3,5 | bedGraphPack stdin ${signal_root}.bedGraph
    bedGraphToBigWig ${signal_root}.bedGraph chromSizes.txt ${signal_root}.bw
    set +x

    echo "* Prepare metadata..."
    qc_hotspot=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_hotspot=`qc_metrics.py -n hotspot -f ${bam_root}_hotspot_out.txt`
        qc_spots=`qc_metrics.py -n singleton -f ${narrowPeak_root}_qc.txt -k "peak count" --keypair "peak count"`
        qc_regions=`qc_metrics.py -n singleton -f ${broadPeak_root}_qc.txt -k "hotspot count" --keypair "hotspot count"`
        qc_hotspot=`echo $qc_hotspot, \"peak_counts\": { $qc_spots, $qc_regions }`
    fi
    # All qc to one file:
    echo "===== hotspot out ====="                > ${bam_root}_hotspot_qc.txt
    cat ${bam_root}_hotspot_out.txt              >> ${bam_root}_hotspot_qc.txt
    echo " "                                     >> ${bam_root}_hotspot_qc.txt
    echo "===== hotspot narrowPeak count ====="  >> ${bam_root}_hotspot_qc.txt
    echo " "                                     >> ${bam_root}_hotspot_qc.txt
    cat ${narrowPeak_root}_qc.txt                >> ${bam_root}_hotspot_qc.txt
    echo "===== hotspot broadPeak count ====="   >> ${bam_root}_hotspot_qc.txt
    cat ${broadPeak_root}_qc.txt                 >> ${bam_root}_hotspot_qc.txt
    
    echo "* Upload results..."
    bed_hotspot_narrowPeak=$(dx upload ${narrowPeak_root}.bed --details "{ $qc_spots }"   --property SW="$versions" --brief)
    bed_hotspot_broadPeak=$(dx upload ${broadPeak_root}.bed   --details "{ $qc_regions }" --property SW="$versions" --brief)
    bb_hotspot_narrowPeak=$(dx upload ${narrowPeak_root}.bb   --details "{ $qc_spots }"   --property SW="$versions" --brief)
    bb_hotspot_broadPeak=$(dx upload ${broadPeak_root}.bb     --details "{ $qc_regions }" --property SW="$versions" --brief)
    bw_hotspot_signal=$(dx upload ${signal_root}.bw           --details "{ $qc_hotspot }" --property SW="$versions" --brief)
    bam_hotspot_qc=$(dx upload ${bam_root}_hotspot_qc.txt --details "{ $qc_hotspot }" --property SW="$versions" --brief)

    dx-jobutil-add-output bed_hotspot_narrowPeak "$bed_hotspot_narrowPeak" --class=file
    dx-jobutil-add-output bed_hotspot_broadPeak "$bed_hotspot_broadPeak" --class=file
    dx-jobutil-add-output bb_hotspot_narrowPeak "$bb_hotspot_narrowPeak" --class=file
    dx-jobutil-add-output bb_hotspot_broadPeak "$bb_hotspot_broadPeak" --class=file
    dx-jobutil-add-output bw_hotspot_signal "$bw_hotspot_signal" --class=file
    dx-jobutil-add-output bam_hotspot_qc "$bam_hotspot_qc" --class=file
    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
