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
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
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
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Fastq${concat} file: '${outfile_name}.fq.gz'"
    reads_root=${outfile_name}
    bam_root="${reads_root}"
    if [ -f /usr/bin/parse_property.py ]; then
        new_root=`parse_property.py --job "${DX_JOB_ID}" --root_name --quiet`
        if [ "$new_root" != "" ]; then
            bam_root="${new_root}"
        fi
    fi

    bwa_ix_root=`dx describe "$bwa_index" --name`
    bwa_ix_root=${bwa_ix_root%.tar.gz}
    bwa_ix_root=${bwa_ix_root%.tgz}
    echo "* Downloading ${bwa_ix_root}.tgz file..."
    dx download "$bwa_index" -o ${bwa_ix_root}.tgz

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_align_bwa_se.sh ${bwa_ix_root}.tgz ${reads_root}.fq.gz $nthreads $bam_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    scripted_root="${bam_root}_se_bwa"
    bam_root="${scripted_root}_techrep"
    # Add DX/encodeD specific _techrep qualifier
    set -x
    mv ${scripted_root}.bam ${bam_root}.bam 
    mv ${scripted_root}_flagstat.txt ${bam_root}_flagstat.txt 
    mv ${scripted_root}_edwBamStats.txt ${bam_root}_edwBamStats.txt 
    set +x
    echo "-- The named results..."
    ls -l ${bam_root}*

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
    bam_bwa=$(dx upload ${bam_root}.bam --details "{ $qc_aligned }" --property SW="$versions" \
                                        --property mapped_reads="$mapped_reads" --property all_reads="$all_reads" \
                                        --property read_length="$read_len" --brief)
    bam_qc=$(dx upload ${bam_root}_qc.txt --details "{ $qc_aligned }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_bwa "$bam_bwa" --class=file
    dx-jobutil-add-output bam_qc "$bam_qc" --class=file

    dx-jobutil-add-output all_reads "$all_reads" --class=string
    dx-jobutil-add-output mapped_reads "$mapped_reads" --class=string
    dx-jobutil-add-output metadata "$qc_aligned" --class=string

    echo "* Finished."
}
