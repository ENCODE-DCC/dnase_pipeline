#!/bin/bash
# dnase-size-bam.sh  Conditionally reduce the bam size to a target, or no larger than a limit of reads

script_name="dnase-size-bam.sh"
script_ver="0.0.1"

main() {
    # Executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of unsized_bam: '$unsized_bam'"
    echo "* Value of target_size: '$target_size'"
    echo "* Value of upper_limit: '$upper_limit'"

    echo "* Download files..."
    bam_root=`dx describe "$unsized_bam" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_root" -o ${bam_root}_unsized.bam
    echo "* unsized_bam file: '${bam_root}_unsized.bam'"

    echo "* bam_sized will be: '${bam_root}_sized.bam'"

    # Skip the contortions of looking up value in details or properties.  Just run edwBamStats!
    echo "* Running edwBamStats on '${bam_root}_unsized.bam'"
    set -x
    edwBamStats ${bam_root}_unsized.bam ${bam_root}_unsized_stats.txt
    set +x
    reads_unsized=`qc_metrics.py -n edwBamStats -f ${bam_root}_unsized_stats.txt -k readCount`
    
    # Figure out the right size
    reads_sized=$reads_unsized
    if [ $reads_unsized -gt $target_size ]; then
        # only if unsized is more than 1.5 times target
        new_target=`expr $target_size \* 3 / 2`  # integers only
        if [ $reads_unsized -gt $new_target ]; then
            reads_sized=$target_size
        fi
    fi
    if [ $upper_limit -gt 0 ] && [ $reads_sized -gt $upper_limit ]; then
        reads_sized=$upper_limit
    fi
    
    if [ $reads_sized -eq $reads_unsized ]; then
        echo "* Current size of '${bam_root}_unsized.bam' will do"
        # Nothing to be done?
        set -x
        mv ${bam_root}_unsized.bam ${bam_root}_sized.bam
        mv ${bam_root}_unsized_stats.txt ${bam_root}_sized.bam_stats.txt
        set +x
        # Output unsized as sized?
    else
        echo "* Sampling $reads_sized reads from '${bam_root}_unsized.bam'"
        set -x
        edwBamStats -sampleBamSize=$reads_sized -u4mSize=$reads_sized -sampleBam=${bam_root}_sized.bam \
                                                              ${bam_root}_unsized.bam ${bam_root}_sized.bam_stats.txt
        set +x
        # What would be nice is to move $unsized_bam to ${bam_root}_full.bam 
        # then upload ${bam_root}_sized.bam to ${bam_root}.bam
        # This would allow blindly using ${bam_root}.bam going forward.
        # Instead, this app will always return a bam which may simply be a copy of the old one.
    fi

    echo "* Prepare metadata for sized bam..."
    qc_sized=''
    final_reads_count=$reads_sized
    reads_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_sized=`qc_metrics.py -n edwBamStats -f ${bam_root}_sized.bam_stats.txt`
        final_reads_count=`qc_metrics.py -n edwBamStats -f ${bam_root}_sized_stats.txt -k readCount`
        reads_len=`qc_metrics.py -n edwBamStats -f ${bam_root}_sized_stats.txt -k readSizeMean`
    fi
    
    echo "* Upload results..."
    bam_sized=$(dx upload ${bam_root}_sized.bam --details "{ $qc_sized }" --property QC="{ $qc_sized }" \
                                                --property reads="$final_reads_count" --property read_length="$reads_len" \
                                                --property SW="$versions" --brief)
    bam_sized_qc=$(dx upload ${bam_root}_sized.bam_stats.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_sized    "$bam_sized"    --class=file
    dx-jobutil-add-output bam_sized_qc "$bam_sized_qc" --class=file
    
    dx-jobutil-add-output final_read_count "$final_reads_count" --class=string
    dx-jobutil-add-output metadata "$versions" --class=string
    
    echo "* Finished."
}
