#!/bin/bash
# fastq-stats.sh

script_name="fastq-stats.sh"
script_ver="0.2.0"

main() {
    # Executable in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of reads: '$reads'"

    #echo "* Download files..."
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!reads[@]}
    do
        filename=`dx describe "${reads[$ix]}" --name | cut -d'.' -f1`
        file_root=${filename%.fastq.gz}
        file_root=${filename%.fq.gz}
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
    
    echo "* Running fastqStatsAndSubsample..."
    set -x
    fastqStatsAndSubsample -smallOk -seed=12345 ${reads_root}.fq.gz ${reads_root}_qc.txt ${reads_root}_sample.fq
    set +x

    echo "* Prepare metadata json..."
    meta=""
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n fastqStatsAndSubsample -f ${reads_root}_qc.txt`
        for ix in ${!reads[@]}
        do
            filename=`dx describe "${reads[$ix]}" --name`
            echo "* Set QC property in ${filename}..."
            set -x
            $(dx set_properties "${reads[$ix]}" QC="{ $meta }")
            set +x
        done
    fi

    echo "* Upload results..."
    fastq_qc=$(dx upload ${reads_root}_qc.txt --details "{ $meta }" --property QC="{ $meta }" --property SW="$versions" --brief)
    #gzip ${reads_root}_sample.fq
    #fastq_sample=$(dx upload ${reads_root}_sample.fq.gz --brief)

    dx-jobutil-add-output fastq_qc "$fastq_qc" --class=file
    #dx-jobutil-add-output fastq_sample "$fastq_sample" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
