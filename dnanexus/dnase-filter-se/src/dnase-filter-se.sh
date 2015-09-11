#!/bin/bash
# dnase-filter-se.sh - Filter bam (single-end) for the ENCODE DNase-seq pipeline.

main() {
    # executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_bwa: '$bam_bwa'"
    echo "* Value of map_thresh: '$map_thresh'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    # expecting *_concat_bwa_biorep.bam
    bam_bwa_root=`dx describe "$bam_bwa" --name`
    #bam_bwa_root=${bam_bwa_root%_concat_bwa_biorep.bam}
    #bam_bwa_root=${bam_bwa_root%_bwa_biorep.bam}
    #bam_bwa_root=${bam_bwa_root%_biorep.bam}
    #bam_bwa_root=${bam_bwa_root%_bwa.bam}
    bam_bwa_root=${bam_bwa_root%.bam}
    dx download "$bam_bwa" -o ${bam_bwa_root}.bam
    echo "* bam_bwa file: '${bam_bwa_root}.bam'"
    bam_filtered_root="${bam_bwa_root}_filtered"
    
    echo "* Filter on threashold..."
    # -F 1804 means not:  0111 0000 1100
    #       4 read unmapped
    #       8 mate unmapped
    #     256 not primary alignment
    #     512 read fails platform/vendor quality checks
    #    1024 read is PCR or optical duplicate
    # -F 780 means:  0011 0000 1100 not: 4,8,256,512
    set -x
    samtools view -F 780 -q ${map_thresh} -u ${bam_bwa_root}.bam | \
            samtools sort -@ $nthreads -m 6G -f - ${bam_filtered_root}.sam
    samtools view -hb ${bam_filtered_root}.sam > ${bam_filtered_root}.bam
    samtools index ${bam_filtered_root}.bam
    set +x

    echo "* Collect filtered bam stats..."
    set -x
    samtools flagstat ${bam_filtered_root}.bam > ${bam_filtered_root}_flagstat.txt
    samtools stats ${bam_filtered_root}.bam > ${bam_filtered_root}_samstats.txt
    head -3 ${bam_filtered_root}_samstats.txt
    grep ^SN ${bam_filtered_root}_samstats.txt | cut -f 2- > ${bam_filtered_root}_samstats_summary.txt
    set +x

    echo "* Prepare metadata for filtered bam..."
    qc_filtered=''
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_filtered=`qc_metrics.py -n samtools_flagstats -f ${bam_filtered_root}_flagstat.txt`
        reads=`qc_metrics.py -n samtools_flagstats -f ${bam_filtered_root}_flagstat.txt -k total`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_filtered_root}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_filtered_root}_samstats_summary.txt -k "average length"`
        qc_filtered=`echo $qc_filtered, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat ====="   > ${bam_filtered_root}_qc.txt
    cat ${bam_filtered_root}_flagstat.txt >> ${bam_filtered_root}_qc.txt
    echo " "                              >> ${bam_filtered_root}_qc.txt
    echo "===== samtools stats ====="     >> ${bam_filtered_root}_qc.txt
    cat ${bam_filtered_root}_samstats.txt >> ${bam_filtered_root}_qc.txt
    
    echo "* Upload results..."
    bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "{ $qc_filtered }" --property SW="$versions" \
                                                      --property reads="$reads" --property read_length="$read_len" --brief)
    bam_filtered_qc=$(dx upload ${bam_filtered_root}_qc.txt --details "{ $qc_filtered }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_filtered "$bam_filtered" --class=file
    dx-jobutil-add-output bam_filtered_qc "$bam_filtered_qc" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
