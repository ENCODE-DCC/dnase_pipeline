#!/bin/bash
# dnase-merge-bams.sh - Merge two or more bams for the ENCODE DNase-seq pipeline.

script_name="dnase-merge-bams.sh"
script_ver="0.2.1"

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of bam_set: '$bam_set'"
    echo "* Value of nthreads: '$nthreads'"

    outfile_name=""
    merged=""
    rm -f concat.fq
    for ix in ${!bam_set[@]}
    do
        filename=`dx describe "${bam_set[$ix]}" --name`
        file_root=${filename%_bwa.bam}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${merged}" == "" ]; then
                outfile_name="${outfile_name}_bwa_merged" 
                merged="s merged as"
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
    
    # TODO: sorting needed?
    echo "* Sorting merged bam..."
    set -x
    samtools sort -@ $nthreads -m 50G -f sofar.bam sorted.bam
    samtools view -hb sorted.bam > sofar.bam 
    set +x
    
    # At this point there is a 'sofar.bam' with one or more input bams
    if [ "${merged}" == "" ]; then
        # Needs to end in '_merged.bam' for pattern matching
        outfile_name="${file_root}_bwa_not_merged"
        mv sofar.bam ${outfile_name}.bam.bam
        echo "* Only one input file, no merging required."
    else
        mv sofar.bam ${outfile_name}.bam
        echo "* Files merged into '${outfile_name}.bam'"
    fi 

    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${outfile_name}.bam > ${outfile_name}_bam_qc.txt
    samtools stats ${outfile_name}.bam > ${outfile_name}_bam_qc_full.txt
    head -3 ${outfile_name}_bam_qc_full.txt
    grep ^SN ${outfile_name}_bam_qc_full.txt | cut -f 2- > ${outfile_name}_bam_qc_summary.txt
    set +x


    echo "* Prepare metadata..."
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n samtools_flagstats -f ${outfile_name}_bam_qc.txt`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${outfile_name}_bam_qc_summary.txt`
        qc_stats=`echo $qc_stats, $meta`
    fi

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_merged=$(dx upload ${outfile_name}.bam --details "{ $qc_stats }" --property QC="{ $qc_stats }" --property SW="$versions" --brief)
    bam_merged_qc=$(dx upload ${outfile_name}_bam_qc.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_merged "$bam_merged" --class=file
    dx-jobutil-add-output bam_merged_qc "$bam_merged_qc" --class=file
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
