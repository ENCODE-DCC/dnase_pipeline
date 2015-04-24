#!/bin/bash
# fastq-stats 0.0.1

main() {
    # Executable in resources/usr/bin
    
    echo "*****"
    echo "* Running: fastq-stats.sh v0.0.1"
    echo "* fastqStatsAndSubsample version: "`fastqStatsAndSubsample 2>&1 | grep "fastqStatsAndSubsample v" | awk '{print $2}'`
    echo "*****"

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
    meta=`echo \"fastqStatsAndSubsample\": { `
    #                  readCount 1000000
    var=`grep -w       readCount ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta \"readCount\": $var`
    #                  baseCount 76000000
    var=`grep -w       baseCount ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"baseCount\": $var`
    #                  sampleCount 100000
    var=`grep -w       sampleCount ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"sampleCount\": $var`
    #                  basesInSample 0
    var=`grep -w       basesInSample ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"basesInSample\": $var`
    #                  readSizeMean 76
    var=`grep -w       readSizeMean ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"readSizeMean\": $var`
    #                  readSizeStd 0
    var=`grep -w       readSizeStd ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"readSizeStd\": $var`
    #                  readSizeMin 76
    var=`grep -w       readSizeMin ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"readSizeMin\": $var`
    #                  readSizeMax 76
    var=`grep -w       readSizeMax ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"readSizeMax\": $var`
    #                  qualMean 32.0269
    var=`grep -w       qualMean ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualMean\": $var`
    #                  qualStd 10.048
    var=`grep -w       qualStd ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualStd\": $var`
    #                  qualMin 2
    var=`grep -w       qualMin ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualMin\": $var`
    #                  qualMax 40
    var=`grep -w       qualMax ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualMax\": $var`
    #                  qualType solexa
    var=`grep -w       qualType ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualType\": \"$var\"`
    #                  qualZero 64
    var=`grep -w       qualZero ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualZero\": $var`
    #                  atRatio 0.535794
    var=`grep -w       atRatio ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"atRatio\": $var`
    #                  aRatio 0.266893
    var=`grep -w       aRatio ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"aRatio\": $var`
    #                  cRatio 0.228095
    var=`grep -w       cRatio ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"cRatio\": $var`
    #                  gRatio 0.235847
    var=`grep -w       gRatio ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"gRatio\": $var`
    #                  tRatio 0.268596
    var=`grep -w       tRatio ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"tRatio\": $var`
    #                  nRatio 0.000568763
    var=`grep -w       nRatio ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"nRatio\": $var`
    #                  posCount 76
    var=`grep -w       posCount ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"posCount\": $var`
    #                  qualPos 35.3192,34.7235,31.5246,34.5956,...,
    var=`grep -w       qualPos ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"qualPos\": [ ${var%,} ]`
    #                  aAtPos 0.2425,0.243239,0.24085,0.282157,...,
    var=`grep -w       aAtPos ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"aAtPos\": [ ${var%,} ]`
    #                  cAtPos 0.222493,0.206261,0.219948,0.20502,...,
    var=`grep -w       cAtPos ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"cAtPos\": [ ${var%,} ]`
    #                  gAtPos 0.324532,0.292566,0.254245,0.278443,...,
    var=`grep -w       gAtPos ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"gAtPos\": [ ${var%,} ]`
    #                  tAtPos 0.208397,0.257934,0.249656,0.234379,...,
    var=`grep -w       tAtPos ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"tAtPos\": [ ${var%,} ]`
    #                  nAtPos 0.002078,0,0.035301,1e-06,0,0,0,0,0,...,
    var=`grep -w       nAtPos ${fastq_root}_qc.txt | awk '{print $2}'`
    meta=`echo $meta,\"nAtPos\": [ ${var%,} ] }`

    echo "* Upload results..."
    fastq_qc=$(dx upload ${fastq_root}_qc.txt --brief)
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
