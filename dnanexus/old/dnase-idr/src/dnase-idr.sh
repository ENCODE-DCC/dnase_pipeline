#!/bin/bash
# dnase-idr.sh

main() {
    echo "* Installing Anaconda3-2.2.0 (python3.4.3, numpy-1.9.2 matplotlib-1.4.3..."
    set -x
    wget https://repo.continuum.io/archive/Anaconda3-2.2.0-Linux-x86_64.sh >> ../install.log 2>&1
    bash Anaconda3-2.2.0-Linux-x86_64.sh -b
    ana_bin=`echo ~/anaconda3/bin`
    # python symlink will interfere with python2.7
    rm ${ana_bin}/python
    #ls ${ana_bin}
    export PATH=${ana_bin}:$PATH
    python3 -V 2>&1 | tee -a install.log
    #sudo python3 -V 2>&1 | tee -a install.log
    #sudo ${ana_bin}/python3 -V 2>&1 | tee -a install.log
    set +x

    echo "* Installing idr..."
    set -x
    wget https://github.com/nboley/idr/archive/2.0.2.tar.gz -O idr.tgz >> ../install.log 2>&1
    mkdir idr
    tar -xzf idr.tgz -C idr --strip-components=1
    cd idr
    # sudo does not see python3 so it requires ana_bin path
    sudo ${ana_bin}/python3 setup.py install >> ../install.log 2>&1
    cd ..
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "Value of peaks_a:     '$peaks_a'"
    echo "Value of peaks_b:     '$peaks_b'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    rank_col="signal.value"
    format="narrowPeak"
    suffix=""
    peaks_a_file=`dx describe "$peaks_a" --name`
    peaks_a_root=${peaks_a_file%.gz}
    peaks_a_root=${peaks_a_root%.bed}
    if [[ "$peaks_a_root" == *_peaks ]]; then
        suffix="peaks"
        rank_col="signal.value"
        format="narrowPeak"
    elif [[ "$peaks_a_root" == *_hotspots ]]; then
        suffix="hotspots"
        rank_col="q.value"
        format="broadPeak"
    fi
    peaks_a_root=${peaks_a_root%_peaks}
    peaks_a_root=${peaks_a_root%_hotspots}
    peaks_a_root=${peaks_a_root%_filtered}
    peaks_a_root=${peaks_a_root%_biorep}
    peaks_a_root=${peaks_a_root%_bwa}
    peaks_a_root=${peaks_a_root%_se}
    peaks_a_root=${peaks_a_root%_pe}
    dx download "$peaks_a" -o $peaks_a_file
    echo "* First bed: '${peaks_a_file}'"

    suffix2=""
    peaks_b_file=`dx describe "$peaks_b" --name`
    peaks_b_root=${peaks_a_file%.gz}
    peaks_b_root=${peaks_b_root%.bed}
    if [[ "$peaks_a_root" == *_peaks ]]; then
        suffix2="peaks"
    elif [[ "$peaks_a_root" == *_hotspots ]]; then
        suffix2="hotspots"
    fi
    peaks_b_root=${peaks_b_root%_peaks}
    peaks_b_root=${peaks_b_root%_hotspots}
    peaks_b_root=${peaks_b_root%_filtered}
    peaks_b_root=${peaks_b_root%_biorep}
    peaks_b_root=${peaks_b_root%_bwa}
    peaks_b_root=${peaks_b_root%_se}
    peaks_b_root=${peaks_b_root%_pe}
    dx download "$peaks_b" -o $peaks_b_file
    echo "* Second bed: '${peaks_b_file}'"

    dx download "$chrom_sizes" -o chrom.sizes

    if [ "$suffix" != "$suffix2" ]; then
        echo "* WARNING: ${peaks_a_file} does not match ${peaks_b_file}: '$suffix' != '$suffix2'."
        suffix="mixed"
    fi
    idr_root="${peaks_a_root}_to_${peaks_b_root}_${suffix}_dnase_idr"
    echo "* IDR root: '"$idr_root"'"

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_idr.sh $peaks_a_file $peaks_b_file $format $rank_col chrom.sizes $idr_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    
    echo "* Conpressing bed file..."
    set -x
    pigz ${idr_root}.bed
    set +x

    echo "* Prepare metadata..."
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n IDR_summary -f idr_summary.txt`
    fi
    
    echo "* Upload results..."
    idr_bed=$(dx upload ${idr_root}.bed.gz --details="{ $qc_stats }" --property SW="$versions" --brief)
    idr_bb=$(dx upload ${idr_root}.bb      --details="{ $qc_stats }" --property SW="$versions" --brief)
    idr_png=$(dx upload ${idr_root}.png    --details="{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output idr_bed "$idr_bed" --class=file
    dx-jobutil-add-output idr_bb "$idr_bb" --class=file
    dx-jobutil-add-output idr_png "$idr_png" --class=file
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
