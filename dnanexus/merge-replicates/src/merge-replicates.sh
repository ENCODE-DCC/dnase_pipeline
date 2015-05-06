#!/bin/bash
# merge-replicates.sh

script_name="merge-replicates.sh"
script_ver="0.2.0"

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
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
    dx download "$bam_A" -o bam_B.bam
    echo "* bam_B file: '${bam_B_root}_filtered.bam'"
    bam_pooled_root="${bam_A_root}_${bam_B_root}_pooled"
    echo "* bam_pooled will be: '${bam_pooled_root}.bam'"

    peaks_A_root=`dx describe "$peaks_A" --name`
    peaks_A_root=${peaks_A_root%_narrowPeak_hotspot.bb}
    dx download "$peaks_A" -o peaks_A.bb
    echo "* peaks_A file: '${peaks_A_root}_narrowPeak_hotspot.bb'"

    peaks_B_root=`dx describe "$peaks_B" --name`
    peaks_B_root=${peaks_B_root%_narrowPeak_hotspot.bb}
    dx download "$peaks_B" -o peaks_B.bb
    echo "* peaks_B file: '${peaks_B_root}_narrowPeak_hotspot.bb'"
    peaks_root="${peaks_A_root}_${peaks_B_root}" 
    peaks_merged_root="${peaks_root}_merged_narrowPeak" 
    echo "* peaks_merged will be: '${peaks_merged_root}.bed/.bb'"

    signal_A_root=`dx describe "$signal_A" --name`
    signal_A_root=${signal_A_root%_signal_hotspot.bw}
    dx download "$signal_A" -o signal_A.bw
    echo "* signal_A file: '${signal_A_root}_signal_hotspot.bw'"

    signal_B_root=`dx describe "$signal_B" --name`
    signal_B_root=${signal_B_root%_signal_hotspot.bw}
    dx download "$signal_B" -o signal_B.bw
    echo "* signal_B file: '${signal_B_root}_signal_hotspot.bw'"
    signal_root="${signal_A_root}_${signal_B_root}_signal"

    dx download "$chrom_sizes" -o chromSizes.txt

    echo "* Merging bams..."
    set -x
    samtools cat bam_A.bam bam_B.bam > ${bam_pooled_root}.bam
    set +x

    echo "* Merging peaks..."
    set -x
    bigBedToBed peaks_A.bb peaks_A.bed
    bigBedToBed peaks_B.bb peaks_B.bed
    cat peaks_A.bed peaks_B.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > ${peaks_merged_root}.bed
    bedToBigBed ${peaks_merged_root}.bed chromSizes.txt ${peaks_merged_root}.bb
    set +x

    echo "* Correlating signals..."
    set -x
    bigWigCorrelate -restrict=${peaks_merged_root}.bb signal_A.bw signal_B.bw > ${signal_root}_corr_qc.txt
    set +x

    echo "* Comparing peak overlaps..."
    set -x
    edwComparePeaks peaks_A.bb peaks_B.bb ${peaks_root}_overlap_qc.txt
    set +x

    # TODO: Collect QC from: ${signal_root}_corr_qc.txt and ${peaks_root}_overlap_qc.txt?
    
    echo "* Prepare metadata..."
    qc_peaks=''  
    qc_signal='' 
    if [ -f /usr/bin/qc_metrics.py ]; then
        #qc_signal=`qc_metrics.py -n bigWigCorrelate -f ${signal_root}_corr_qc.txt`
        qc_signal=`qc_metrics.py -n singleton -f ${signal_root}_corr_qc.txt -k "bigWigCorrelate" --keypair "bigWigCorrelate"`
        qc_peaks=`qc_metrics.py -n edwComparePeaks -f ${peaks_root}_overlap_qc.txt`
        qc_signal=`echo $qc_signal, $qc_peaks`
    fi

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_pooled=$(dx upload ${bam_pooled_root}.bam   --details "{ $qc_signal }" --property QC="{ $qc_signal }" --property SW="$versions" --brief)
    bed_merged=$(dx upload ${peaks_merged_root}.bed --details "{ $qc_peaks }" --property QC="{ $qc_peaks }" --property SW="$versions" --brief)
    bb_merged=$(dx upload ${peaks_merged_root}.bb   --details "{ $qc_peaks }" --property QC="{ $qc_peaks }" --property SW="$versions" --brief)
    signal_corr_qc=$(dx upload ${signal_root}_corr_qc.txt --property SW="$versions" --brief)
    peaks_overlap_qc=$(dx upload ${peaks_root}_overlap_qc.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_pooled "$bam_pooled" --class=file
    dx-jobutil-add-output bed_merged "$bed_merged" --class=file
    dx-jobutil-add-output bb_merged "$bb_merged" --class=file
    dx-jobutil-add-output signal_corr_qc "$signal_corr_qc" --class=file
    dx-jobutil-add-output peaks_overlap_qc "$peaks_overlap_qc" --class=file
    dx-jobutil-add-output metadata "{ $qc_signal }" --class=string

    echo "* Finished."
}
