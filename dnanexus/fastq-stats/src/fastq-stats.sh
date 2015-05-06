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

    echo "* Value of fastq_file: '$fastq_file'"

    #echo "* Download files..."
    if [ ${#fastq_file} -gt 1 ]; then
        outfile_name="concat"
        rm -f concat.fq
        for ix in ${!fastq_file[@]}
        do
            filename=`dx describe "${fastq_file[$ix]}" --name | cut -d'.' -f1`
            file_root=${filename%.fastq.gz}
            file_root=${filename%.fq.gz}
            outfile_name="${file_root}_${outfile_name}"
            echo "* Downloading and concatenating ${file_root}.fq.gz file..."
            dx download "${fastq_file[$ix]}" -o - | gunzip >> concat.fq
        done
        mv concat.fq ${outfile_name}.fq
        echo "* Gzipping file..."
        gzip ${outfile_name}.fq
        echo "* Fastqs concatenated as: '${outfile_name}.fq.gz'"
        fastq_root=${outfile_name}
    else
        echo "* Download files..."
        fastq_root=`dx describe "$fastq_file" --name`
        fastq_root=${fastq_root%.fastq.gz}
        fastq_root=${fastq_root%.fq.gz}
        dx download "$fastq_file" -o ${fastq_root}.fq.gz
        echo "* Fastq file: '${fastq_root}.fq.gz'"
    fi

    echo "* Running fastqStatsAndSubsample..."
    set -x
    fastqStatsAndSubsample -smallOk -seed=12345 ${fastq_root}.fq.gz ${fastq_root}_qc.txt ${fastq_root}_sample.fq
    set +x

    echo "* Prepare metadata json..."
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n fastqStatsAndSubsample -f ${fastq_root}_qc.txt`
    fi

    echo "* Upload results..."
    fastq_qc=$(dx upload ${fastq_root}_qc.txt --details "{ $meta }" --property QC="{ $meta }" --property SW="$versions" --brief)
    #gzip ${fastq_root}_sample.fq
    #fastq_sample=$(dx upload ${fastq_root}_sample.fq.gz --brief)

    dx-jobutil-add-output fastq_qc "$fastq_qc" --class=file
    #dx-jobutil-add-output fastq_sample "$fastq_sample" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    if [ ${#fastq_file} eq 1 ]; then
        dx set-properties "${fastq_file}" QC="$meta"
    fi

    echo "* Finished."
}
