#!/bin/bash
# merge-replicates.sh

script_name="merge-replicates.sh"
script_ver="0.0.1"

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi
    #echo "*****"
    #echo "* Running: merge-replicates.sh v0.0.1"
    #echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    #echo "* bedtools version: "`bedtools --version 2>&1 | awk '{print $2}'`
    #echo "* bigBedToBed version: "`bigBedToBed 2>&1 | grep "bigBedToBed v" | awk '{print $2}'`
    #echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{print $2$3}'`
    #echo "* bigWigCorrelate version: (unversioned) "`bigWigCorrelate 2>&1 | grep "bigWigCorrelate -" | awk '{print $3,$4,$5}'`
    #echo "* edwComparePeaks version: (unversioned) "`edwComparePeaks 2>&1 | grep "edwComparePeaks -" | awk '{print $3,$4,$5,$6}'`
    #echo "*****"

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
    echo "* peaks_merged will be: '${peaks_merged_root}_merged_narrowPeak.bed/.bb'"

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
    cat peaks_A.bed peaks_B.bed | sort -k1,1 -k2,2n - | bedtools merge -i - > ${peaks_root}_merged_narrowPeak.bed
    bedToBigBed ${peaks_root}_merged_narrowPeak.bed chromSizes.txt ${peaks_root}_merged_narrowPeak.bb
    set +x

    echo "* Correlating signals..."
    set -x
    bigWigCorrelate -restrict=${peaks_root}_merged_narrowPeak.bb signal_A.bw signal_B.bw 1> ${signal_root}_corr_qc.txt
    set +x

    echo "* Comparing peak overlaps..."
    set -x
    edwComparePeaks peaks_A.bb peaks_B.bb ${peaks_root}_overlap_qc.txt
    set +x

    # TODO: Collect QC from: ${signal_root}_corr_qc.txt and ${peaks_root}_overlap_qc.txt?
    
    echo "* Prepare metadata..."
    meta=`echo \"bigWigCorrelate\": { `
    ## # 2142994 + 0 in total (QC-passed reads + QC-failed reads)
    ## var=`grep "QC-passed reads" ${bam_filtered_root}_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    ## meta=`echo $meta $var`
    ## # 0 + 0 duplicates
    ## var=`grep -w duplicates ${bam_filtered_root}_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    ## meta=`echo $meta, $var`
    ## # 2046212 + 0 mapped (95.48%:-nan%)
    ## var=`grep -w mapped ${bam_filtered_root}_qc.txt | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    ## meta=`echo $meta, $var`
    ## var=`grep -w mapped ${bam_filtered_root}_qc.txt | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    ## meta=`echo $meta, \"mapped_pct\": \"$var\"`
    meta=`echo $meta }`

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    #bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "{ $meta }" --property QC="{ $meta }" --brief)
    bam_pooled=$(dx upload ${bam_pooled_root}.bam --brief)
    bed_merged=$(dx upload ${peaks_root}_merged_narrowPeak.bed --brief)
    bb_merged=$(dx upload ${peaks_root}_merged_narrowPeak.bb --brief)
    signal_corr_qc=$(dx upload ${signal_root}_corr_qc.txt --brief)
    peaks_overlap_qc=$(dx upload ${peaks_root}_overlap_qc.txt --brief)

    dx-jobutil-add-output bam_pooled "$bam_pooled" --class=file
    dx-jobutil-add-output bed_merged "$bed_merged" --class=file
    dx-jobutil-add-output bb_merged "$bb_merged" --class=file
    dx-jobutil-add-output signal_corr_qc "$signal_corr_qc" --class=file
    dx-jobutil-add-output peaks_overlap_qc "$peaks_overlap_qc" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
