#!/bin/bash
# dnase-merge-bams.sh - Merge two or more bams for the ENCODE DNase-seq pipeline.

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_set: '$bam_set'"
    echo "* Value of nthreads: '$nthreads'"
    
    outfile_name=""
    merged=""
    tech_reps=""
    pe_se=""
    rm -f concat.fq
    for ix in ${!bam_set[@]}
    do
        file_root=`dx describe "${bam_set[$ix]}" --name`
        file_root=${file_root%_bwa_techrep.bam}
        file_root=${file_root%_bwa.bam}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${merged}" == "" ]; then
                outfile_name="${outfile_name}_bwa_biorep" 
                merged="s merged as"
            fi
        fi
        if [ -f /usr/bin/parse_property.py ]; then
            # If even one is se then all are se
            sans_se=${file_root%_se}
            if [ "${sans_se}" != "${file_root}" ]; then
                pe_se="_se"
            elif [ "${pe_se}" == "" ]; then
                # If it can be determined to be pe then default it to that
                sans_pe=${file_root%_pe}
                if [ "${sans_pe}" != "${file_root}" ]; then
                    pe_se="_pe"
                fi
            fi
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
        outfile_name="${exp_id}_${tech_reps}${pe_se}_bwa_biorep"
    fi
    echo "* Merged alignments file will be: '${outfile_name}.bam'"
    
    # TODO: sorting needed?
    echo "* Sorting merged bam..."
    set -x
    samtools sort -@ $nthreads -m 6G -f sofar.bam sorted.bam
    samtools view -hb sorted.bam > sofar.bam 
    set +x
    
    # At this point there is a 'sofar.bam' with one or more input bams
    if [ "${merged}" == "" ]; then
        outfile_name="${file_root}_bwa_biorep"
        mv sofar.bam ${outfile_name}.bam
        echo "* Only one input file, no merging required."
    else
        mv sofar.bam ${outfile_name}.bam
        echo "* Files merged into '${outfile_name}.bam'"
    fi 

    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${outfile_name}.bam > ${outfile_name}_flagstat.txt
    samtools stats ${outfile_name}.bam > ${outfile_name}_samstats.txt
    head -3 ${outfile_name}_samstats.txt
    grep ^SN ${outfile_name}_samstats.txt | cut -f 2- > ${outfile_name}_samstats_summary.txt
    set +x


    echo "* Prepare metadata..."
    qc_stats=''
    reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n samtools_flagstats -f ${outfile_name}_flagstat.txt`
        reads=`qc_metrics.py -n samtools_flagstats -f ${outfile_name}_flagstat.txt -k mapped`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${outfile_name}_samstats_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${outfile_name}_samstats_summary.txt -k "average length"`
        qc_stats=`echo $qc_stats, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat =====" > ${outfile_name}_qc.txt
    cat ${outfile_name}_flagstat.txt    >> ${outfile_name}_qc.txt
    echo " "                            >> ${outfile_name}_qc.txt
    echo "===== samtools stats ====="   >> ${outfile_name}_qc.txt
    cat ${outfile_name}_samstats.txt    >> ${outfile_name}_qc.txt

    echo "* Upload results..."
    bam_biorep=$(dx upload ${outfile_name}.bam --details "{ $qc_stats }" --property SW="$versions" \
                                               --property reads="$reads" --property read_length="$read_len" --brief)
    bam_biorep_qc=$(dx upload ${outfile_name}_qc.txt --details "{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_biorep "$bam_biorep" --class=file
    dx-jobutil-add-output bam_biorep_qc "$bam_biorep_qc" --class=file

    dx-jobutil-add-output reads "$reads" --class=string
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
