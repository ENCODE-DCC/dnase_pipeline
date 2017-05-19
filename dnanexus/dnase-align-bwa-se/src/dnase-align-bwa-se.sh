#!/bin/bash
# dnase-align-bwa-se.sh

main() {
    # Executable in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of reads: '$reads'"
    echo "* Value of bwa_index: '$bwa_index'"
    echo "* Value of barcode:   '$barcode'"
    #echo "* Value of UMI:       '$umi'"   # No UMI handling of single-end fastqs... no UMI single-end fastqs expected
    echo "* Value of trim_len:  '$trim_len'"
    echo "* Value of nthreads:  '$nthreads'"

    umi="no"
    if [[ $barcode == SSLIB* ]]; then
        echo "* WARNING: barcode '$barcode' suggests UMI, but single-end is expected to be non-UMI and is being treated as such."
    fi

    #echo "* Download files..."
    exp_rep_root=""
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py -f "'${reads[0]}'" --project "${DX_PROJECT_CONTEXT_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            exp_rep_root="${new_root}"
        fi
    fi
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads[@]}
    do
        file_root=`dx describe "${reads[$ix]}" --name`
        file_root=${file_root%.fastq.gz}
        file_root=${file_root%.fq.gz}
        if [ "${outfile_name}" == "" ]; then
            outfile_name="${file_root}"
        else
            outfile_name="${file_root}_${outfile_name}"
            if [ "${concat}" == "" ]; then
                outfile_name="${outfile_name}_concat" 
                concat="s concatenated as"
            fi
        fi
        echo "* Downloading and concatenating ${file_root}.fq.gz file..."
        dx download "${reads[$ix]}" -o - | gunzip >> concat.fq
    done
    if [ "${concat}" != "" ]; then
        if [ "${exp_rep_root}" != "" ]; then
            outfile_name="${exp_rep_root}"
        elif [ ${#outfile_name} -gt 200 ]; then
            outfile_name="concatenated_reads"
        fi
    fi
    if [ $trim_len -gt 0 ]; then
        echo "* Trimming fastq reads to ${trim_len} bases..."
        outfile_name="${outfile_name}_${trim_len}b" 
        cutadapt -l $trim_len concat.fq > ${outfile_name}.fq
    else
        mv concat.fq ${outfile_name}.fq
    fi
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Fastq${concat} file: '${outfile_name}.fq.gz'"
    reads_root=${outfile_name}
    bam_root="${reads_root}"
    if [ "${exp_rep_root}" != "" ]; then
        bam_root="${exp_rep_root}"
    fi

    bwa_ix_root=`dx describe "$bwa_index" --name`
    bwa_ix_root=${bwa_ix_root%.tar.gz}
    bwa_ix_root=${bwa_ix_root%.tgz}
    echo "* Downloading ${bwa_ix_root}.tgz file..."
    dx download "$bwa_index" -o ${bwa_ix_root}.tgz

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    bam_root="${bam_root}_se_bwa_techrep"
    set -x
    dnase_align_bwa_se.sh ${bwa_ix_root}.tgz ${reads_root}.fq.gz $nthreads $bam_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Prepare metadata json..."
    qc_aligned=''
    all_reads=0
    mapped_reads=0
    read_len=0
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_aligned=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt`
        mapped_reads=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_flagstat.txt -k mapped`
        all_reads=`qc_metrics.py -n edwBamStats -f ${bam_root}_edwBamStats.txt -k readCount`
        read_len=`qc_metrics.py -n edwBamStats -f ${bam_root}_edwBamStats.txt -k readSizeMean`
        meta=`qc_metrics.py -n edwBamStats -f ${bam_root}_edwBamStats.txt`
        qc_aligned=`echo $qc_aligned, $meta`
    fi
    # All qc to one file per target file:
    echo "===== samtools flagstat =====" > ${bam_root}_qc.txt
    cat ${bam_root}_flagstat.txt        >> ${bam_root}_qc.txt
    echo " "                            >> ${bam_root}_qc.txt
    echo "===== edwBamStats ====="      >> ${bam_root}_qc.txt
    cat ${bam_root}_edwBamStats.txt     >> ${bam_root}_qc.txt

    echo "* Upload results..."
    bam_bwa=$(dx upload ${bam_root}.bam --details "{ $qc_aligned }" --property SW="$versions" --property pe_or_se="se" \
                                --property mapped_reads="$mapped_reads" --property all_reads="$all_reads" \
                                --property read_length="$read_len" --property UMI="$umi" --property barcode="$barcode" --brief)
    bam_qc=$(dx upload ${bam_root}_qc.txt --details "{ $qc_aligned }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_bwa "$bam_bwa" --class=file
    dx-jobutil-add-output bam_qc "$bam_qc" --class=file

    dx-jobutil-add-output all_reads "$all_reads" --class=string
    dx-jobutil-add-output mapped_reads "$mapped_reads" --class=string
    dx-jobutil-add-output metadata "$qc_aligned" --class=string

    echo "* Finished."
}
