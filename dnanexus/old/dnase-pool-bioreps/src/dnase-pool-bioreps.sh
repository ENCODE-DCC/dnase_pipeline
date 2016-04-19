#!/bin/bash
# dnase-pool-bioreps.sh - Pools two biological replicates for the ENCODE DNase-seq pipeline.

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_A: '$bam_A'"
    echo "* Value of bam_B: '$bam_B'"
    echo "* Value of peaks_A: '$peaks_A'"
    echo "* Value of peaks_B: '$peaks_B'"
    echo "* Value of signal_A: '$signal_A'"
    echo "* Value of signal_B: '$signal_B'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    bam_A_root=`dx describe "$bam_A" --name`
    bam_A_root=${bam_A_root%_filtered.bam}
    dx download "$bam_A" -o bam_A.bam
    echo "* bam_A file: '${bam_A_root}_filtered.bam'"

    bam_B_root=`dx describe "$bam_B" --name`
    bam_B_root=${bam_B_root%_filtered.bam}
    dx download "$bam_B" -o bam_B.bam
    echo "* bam_B file: '${bam_B_root}_filtered.bam'"
    out_root="${bam_A_root}_${bam_B_root}"
    ### bam_pooled_root="${bam_A_root}_${bam_B_root}_pooled"
    ### echo "* bam_pooled will be: '${bam_pooled_root}.bam'"

    peaks_A_root=`dx describe "$peaks_A" --name`
    peaks_A_root=${peaks_A_root%_narrowPeak_hotspot.bb}
    dx download "$peaks_A" -o peaks_A.bb
    echo "* peaks_A file: '${peaks_A_root}_narrowPeak_hotspot.bb'"

    peaks_B_root=`dx describe "$peaks_B" --name`
    peaks_B_root=${peaks_B_root%_narrowPeak_hotspot.bb}
    dx download "$peaks_B" -o peaks_B.bb
    echo "* peaks_B file: '${peaks_B_root}_narrowPeak_hotspot.bb'"
    ### peaks_root="${peaks_A_root}_${peaks_B_root}" 
    ### peaks_merged_root="${peaks_root}_merged_narrowPeak" 
    ### echo "* peaks_merged will be: '${peaks_merged_root}.bed/.bb'"

    signal_A_root=`dx describe "$signal_A" --name`
    signal_A_root=${signal_A_root%_signal_hotspot.bw}
    dx download "$signal_A" -o signal_A.bw
    echo "* signal_A file: '${signal_A_root}_signal_hotspot.bw'"

    signal_B_root=`dx describe "$signal_B" --name`
    signal_B_root=${signal_B_root%_signal_hotspot.bw}
    dx download "$signal_B" -o signal_B.bw
    echo "* signal_B file: '${signal_B_root}_signal_hotspot.bw'"
    ### signal_root="${signal_A_root}_${signal_B_root}_signal"

    dx download "$chrom_sizes" -o chrom.sizes

    read_len_A=`parse_property.py -f "$bam_A" -p "read_length" --quiet`
    read_len_B=`parse_property.py -f "$bam_B" -p "read_length" --quiet`
    if [ "$read_len_A" == "" ]; then
        echo "* Running edwBamStats on 'bam_A.bam'"
        set -x
        edwBamStats bam_A.bam bam_A_edwBamStats.txt
        set +x
        read_len_A=`qc_metrics.py -n edwBamStats -f bam_A_edwBamStats.txt -k readSizeMean`
    fi
    if [ "$read_len_B" == "" ]; then
        echo "* Running edwBamStats on 'bam_B.bam'"
        set -x
        edwBamStats bam_B.bam bam_B_edwBamStats.txt
        set +x
        read_len_B=`qc_metrics.py -n edwBamStats -f bam_B_edwBamStats.txt -k readSizeMean`
    fi
    if [ "$read_len_A" != "" ] && [ "$read_len_B" != "" ] &&  [ "$read_len_A" != "$read_len_A" ]; then
        echo "* WARNING: Read lengths of two bam files do not match."
    fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_pooled_reps.sh bam_A.bam bam_B.bam peaks_A.bb peaks_B.bb signal_A.bw signal_B.bw chrom.sizes $out_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    bam_pooled_root="${out_root}_pooled"
    peaks_merged_root="${out_root}_merged_narrowPeak" 
    signal_root_root="${out_root}_signal" 
    
    echo "* Compressing bed files..."
    set -x
    gzip ${peaks_merged_root}.bed
    set +x


    echo "* Prepare metadata..."
    qc_pooled=''  
    qc_peaks=''  
    qc_signal='' 
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_pooled=`qc_metrics.py -n edwBamStats -f ${bam_pooled_root}_edwBamStats.txt`
        reads=`qc_metrics.py -n edwBamStats -f ${bam_pooled_root}_edwBamStats.txt -k readCount`
        read_len=`qc_metrics.py -n edwBamStats -f ${bam_pooled_root}_edwBamStats.txt -k readSizeMean`
        #qc_signal=`qc_metrics.py -n bigWigCorrelate -f ${signal_root}_corr_qc.txt`
        qc_signal=`qc_metrics.py -n singleton -f ${signal_root}_corr_qc.txt -k "bigWigCorrelate" --keypair "bigWigCorrelate"`
        qc_peaks=`qc_metrics.py -n edwComparePeaks -f ${out_root}_peaks_overlap_qc.txt`
        qc_pooled=`echo $qc_pooled, $qc_signal, $qc_peaks`
    fi
    # All qc to one file per target file:
    echo "===== edwBamStats ====="          > ${bam_pooled_root}_qc.txt
    cat ${bam_pooled_root}_edwBamStats.txt >> ${bam_pooled_root}_qc.txt
    echo " "                               >> ${bam_pooled_root}_qc.txt
    echo "===== bigWigCorrelate ====="     >> ${bam_pooled_root}_qc.txt
    cat ${signal_root}_corr_qc.txt         >> ${bam_pooled_root}_qc.txt
    echo " "                               >> ${bam_pooled_root}_qc.txt
    echo "===== edwComparePeaks ====="     >> ${bam_pooled_root}_qc.txt
    cat ${out_root}_peaks_overlap_qc.txt   >> ${bam_pooled_root}_qc.txt

    echo "* Upload results..."
    bam_pooled=$(dx upload ${bam_pooled_root}.bam      --details "{ $qc_pooled }" --property SW="$versions" \
                                                       --property reads="$reads" --property read_length="$read_len" --brief)
    bed_merged=$(dx upload ${peaks_merged_root}.bed.gz --details "{ $qc_peaks }" --property SW="$versions" --brief)
    bb_merged=$(dx upload ${peaks_merged_root}.bb      --details "{ $qc_peaks }" --property SW="$versions" --brief)
    pooled_qc=$(dx upload ${bam_pooled_root}_qc.txt    --details "{ $qc_pooled }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_pooled "$bam_pooled" --class=file
    dx-jobutil-add-output bed_merged "$bed_merged" --class=file
    dx-jobutil-add-output bb_merged "$bb_merged" --class=file
    dx-jobutil-add-output pooled_qc "$pooled_qc" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_signal }" --class=string

    echo "* Finished."
}
