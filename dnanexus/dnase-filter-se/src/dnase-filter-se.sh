#!/bin/bash
# dnase-filter-se.sh - Filter bam (single-end) for the ENCODE DNase-seq pipeline.

script_name="dnase-filter-se.sh"
script_ver="0.2.1"

main() {
    # executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of bam_bwa: '$bam_bwa'"
    echo "* Value of map_thresh: '$map_thresh'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    # expecting *_concat_bwa_merged.bam
    bam_bwa_root=`dx describe "$bam_bwa" --name`
    #bam_bwa_root=${bam_bwa_root%_concat_bwa_merged.bam}
    #bam_bwa_root=${bam_bwa_root%_bwa_merged.bam}
    #bam_bwa_root=${bam_bwa_root%_merged.bam}
    #bam_bwa_root=${bam_bwa_root%_bwa.bam}
    bam_bwa_root=${bam_bwa_root%.bam}
    dx download "$bam_bwa" -o ${bam_bwa_root}.bam
    echo "* bam_bwa file: '${bam_bwa_root}.bam'"
    bam_filtered_root="${bam_bwa_root}_filtered"
    
    echo "* Filter on threashold..."
    set -x
    samtools view -F 1804 -q ${map_thresh} -u ${bam_bwa_root}.bam | \
            samtools sort -@ $nthreads -m 5G -f - ${bam_filtered_root}.sam
    samtools view -hb ${bam_filtered_root}.sam > ${bam_filtered_root}.bam
    samtools index ${bam_filtered_root}.bam
    set +x

    echo "* Collect filtered bam stats..."
    set -x
    samtools flagstat ${bam_filtered_root}.bam > ${bam_filtered_root}_qc.txt
    samtools stats ${bam_filtered_root}.bam > ${bam_filtered_root}_qc_full.txt
    head -3 ${bam_filtered_root}_qc_full.txt
    grep ^SN ${bam_filtered_root}_qc_full.txt | cut -f 2- > ${bam_filtered_root}_qc_summary.txt
    set +x

    echo "* Prepare metadata for filtered bam..."
    qc_filtered=''
    reads_filtered=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_filtered=`qc_metrics.py -n samtools_flagstats -f ${bam_filtered_root}_qc.txt`
        reads_filtered=`qc_metrics.py -n samtools_flagstats -f ${bam_filtered_root}_qc.txt -k total`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_filtered_root}_qc_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_filtered_root}_qc_summary.txt --keypair "average length"`
        qc_filtered=`echo $qc_filtered, $meta`
    fi
    
    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "{ $qc_filtered }" --property QC="{ $qc_filtered }" \
                                                      --property reads="$reads_filtered" --property read_length="$reads_len" \
                                                      --property SW="$versions" --brief)
    bam_filtered_qc=$(dx upload ${bam_filtered_root}_qc.txt           --property SW="$versions" --brief)
    bam_filtered_qc_full=$(dx upload ${bam_filtered_root}_qc_full.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_filtered "$bam_filtered" --class=file
    dx-jobutil-add-output bam_filtered_qc "$bam_filtered_qc" --class=file
    dx-jobutil-add-output bam_filtered_qc_full "$bam_filtered_qc_full" --class=file

    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
