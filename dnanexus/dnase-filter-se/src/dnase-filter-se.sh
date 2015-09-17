#!/bin/bash
# dnase-filter-se.sh - Merge and filter bams (single-end) for the ENCODE DNase-seq pipeline.

main() {
    # executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_set: '$bam_set'"
    echo "* Value of map_thresh: '$map_thresh'"
    echo "* Value of nthreads: '$nthreads'"

    merged_bam_root=""
    merged=""
    tech_reps=""
    rm -f concat.fq
    for ix in ${!bam_set[@]}
    do
        file_root=`dx describe "${bam_set[$ix]}" --name`
        file_root=${file_root%_bwa_techrep.bam}
        file_root=${file_root%_bwa.bam}
        if [ "${merged_bam_root}" == "" ]; then
            merged_bam_root="${file_root}"
        else
            merged_bam_root="${file_root}_${merged_bam_root}"
            if [ "${merged}" == "" ]; then
                merged_bam_root="${merged_bam_root}_bwa_biorep" 
                merged="s merged as"
            fi
        fi
        if [ -f /usr/bin/parse_property.py ]; then
            if [ "$exp_id" == "" ]; then
                exp_id=`parse_property.py -f "'${bam_set[$ix]}'" --project "${DX_PROJECT_CONTEXT_ID}" --exp_id`
            fi
            rep_tech=`parse_property.py -f "'${bam_set[$ix]}'" --project "${DX_PROJECT_CONTEXT_ID}" --rep_tech`
            if [ "$rep_tech" != "" ]; then
                if  [ "$tech_reps" != "" ]; then
                    tech_reps="${tech_reps}_${rep_tech}"
                else
                    tech_reps="${rep_tech}"
                fi
            fi
        fi
        echo "* Downloading ${file_root}_bwa.bam file..."
        dx download "${bam_set[$ix]}" -o ${file_root}_bwa.bam
        if [ ! -e sofar.bam ]; then
            mv ${file_root}_bwa.bam sofar.bam
        else
            echo "* Merging..."
            set -x
            samtools cat sofar.bam ${file_root}_bwa.bam > merging.bam
            mv merging.bam sofar.bam
            set +x
        fi
    done
    if [ "$exp_id" != "" ] && [ "$tech_reps" != "" ]; then
        merged_bam_root="${exp_id}_${tech_reps}_se_bwa_biorep"
    fi
    filtered_bam_root="${merged_bam_root}_filtered"
    echo "* Merged alignments file will be: '${merged_bam_root}.bam'"
    echo "* Filtered alignments file will be: '${filtered_bam_root}.bam'"
    
    # At this point there is a 'sofar.bam' with one or more input bams
    if [ "${merged}" == "" ]; then
        merged_bam_root="${file_root}_bwa_biorep"
        mv sofar.bam ${merged_bam_root}.bam
        echo "* Only one input file, no merging required."
    else
        mv sofar.bam ${merged_bam_root}.bam
        echo "* Files merged into '${merged_bam_root}.bam'"
    fi 

    echo "* Filter on threashold..."
    # -F 1804 means not:  0111 0000 1100
    #       4 read unmapped
    #       8 mate unmapped
    #     256 not primary alignment
    #     512 read fails platform/vendor quality checks
    #    1024 read is PCR or optical duplicate
    # -F 780 means:  0011 0000 1100 not: 4,8,256,512
    set -x
    samtools view -F 780 -q ${map_thresh} -u ${merged_bam_root}.bam | \
            samtools sort -@ $nthreads -m 6G -f - ${filtered_bam_root}.sam
    samtools view -hb ${filtered_bam_root}.sam > ${filtered_bam_root}.bam
    samtools index ${filtered_bam_root}.bam
    set +x

    echo "* Collect filtered bam stats..."
    set -x
    samtools flagstat ${filtered_bam_root}.bam > ${filtered_bam_root}_flagstat.txt
    samtools stats ${filtered_bam_root}.bam > ${filtered_bam_root}_samstats.txt
    head -3 ${filtered_bam_root}_samstats.txt
    grep ^SN ${filtered_bam_root}_samstats.txt | cut -f 2- > ${filtered_bam_root}_samstats_summary.txt
    set +x

    echo "* Prepare metadata for filtered bam..."
    qc_filtered=''
    prefiltered_all_reads=0
    prefiltered_mapped_reads=0
    #filtered_all_reads=0
    filtered_mapped_reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_filtered=`qc_metrics.py -n samtools_flagstats -f ${filtered_bam_root}_flagstat.txt`
        prefiltered_all_reads=`qc_metrics.py -n samtools_flagstats -f ${merged_bam_root}_flagstat.txt -k total`
        prefiltered_mapped_reads=`qc_metrics.py -n samtools_flagstats -f ${merged_bam_root}_flagstat.txt -k mapped`
        #filtered_all_reads=`qc_metrics.py -n samtools_flagstats -f ${filtered_bam_root}_flagstat.txt -k total`
        filtered_mapped_reads=`qc_metrics.py -n samtools_flagstats -f ${filtered_bam_root}_flagstat.txt -k mapped`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${filtered_bam_root}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${filtered_bam_root}_samstats_summary.txt -k "average length"`
        qc_filtered=`echo $qc_filtered, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat ====="   > ${filtered_bam_root}_qc.txt
    cat ${filtered_bam_root}_flagstat.txt >> ${filtered_bam_root}_qc.txt
    echo " "                              >> ${filtered_bam_root}_qc.txt
    echo "===== samtools stats ====="     >> ${filtered_bam_root}_qc.txt
    cat ${filtered_bam_root}_samstats.txt >> ${filtered_bam_root}_qc.txt
    
    echo "* Upload results..."
    bam_filtered=$(dx upload ${filtered_bam_root}.bam --details "{ $qc_filtered }" --property SW="$versions" \
                                                      --property prefiltered_all_reads="$prefiltered_all_reads" \
                                                      --property prefiltered_mapped_reads="$prefiltered_mapped_reads" \
                                                      --property filtered_mapped_reads="$filtered_mapped_reads" \
                                                     --property read_length="$read_len" --brief)
     bam_filtered_qc=$(dx upload ${filtered_bam_root}_qc.txt --details "{ $qc_filtered }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_filtered "$bam_filtered" --class=file
    dx-jobutil-add-output bam_filtered_qc "$bam_filtered_qc" --class=file

    dx-jobutil-add-output prefiltered_all_reads "$prefiltered_all_reads" --class=string
    dx-jobutil-add-output prefiltered_mapped_reads "$prefiltered_mapped_reads" --class=string
    dx-jobutil-add-output filtered_mapped_reads "$filtered_mapped_reads" --class=string
    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
