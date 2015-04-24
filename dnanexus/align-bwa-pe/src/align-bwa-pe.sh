#!/bin/bash
# align-bwa-pe 0.1.0

main() {
    # Executable in resources/usr/bin
    
    echo "*****"
    echo "* Running: align-bwa-pe.sh v0.1.0"
    echo "* bwa version: "`bwa 2>&1 | grep Version | awk '{print $2}'`
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "* Value of read1_fq: '$read1_fq'"
    echo "* Value of read2_fq: '$read2_fq'"
    echo "* Value of bwa_index: '$bwa_index'"
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
    if [ "${#read1_fq}" -gt "1" ]; then
        outfile_name="concat"
        rm -f concat.fq
        for ix in ${!read1_fq[@]}
        do
            filename=`dx describe "${read1_fq[$ix]}" --name | cut -d'.' -f1`
            file_root=${filename%.fastq.gz}
            file_root=${filename%.fq.gz}
            outfile_name="${file_root}_${outfile_name}"
            echo "* Downloading and concatenating ${file_root}.fq.gz file..."
            dx download "${read1_fq[$ix]}" -o - | gunzip >> concat.fq
        done
        mv concat.fq ${outfile_name}.fq
        echo "* Gzipping file..."
        gzip ${outfile_name}.fq
        echo "* Fastqs concatenated as: '${outfile_name}.fq.gz'"
        read1_root=${outfile_name}
    else
        read1_root=`dx describe "$read1_fq" --name`
        read1_root=${read1_root%.fastq.gz}
        read1_root=${read1_root%.fq.gz}
        echo "* Download ${read1_root}.fq.gz file..."
        dx download "$read1_fq" -o ${read1_root}.fq.gz
        echo "* Fastq file: '${read1_root}.fq.gz'"
    fi

    if [ "${#read2_fq}" -gt "1" ]; then
        outfile_name="concat"
        rm -f concat.fq
        for ix in ${!read2_fq[@]}
        do
            filename=`dx describe "${read2_fq[$ix]}" --name | cut -d'.' -f1`
            file_root=${filename%.fastq.gz}
            file_root=${filename%.fq.gz}
            outfile_name="${file_root}_${outfile_name}"
            echo "* Downloading and concatenating ${file_root}.fq.gz file..."
            dx download "${read2_fq[$ix]}" -o - | gunzip >> concat.fq
        done
        mv concat.fq ${outfile_name}.fq
        echo "* Gzipping file..."
        gzip ${outfile_name}.fq
        echo "* Fastqs concatenated as: '${outfile_name}.fq.gz'"
        read2_root=${outfile_name}
    else
        read2_root=`dx describe "$read2_fq" --name`
        read2_root=${read2_root%.fastq.gz}
        read2_root=${read2_root%.fq.gz}
        echo "* Downloading ${read2_root}.fq.gz file..."
        dx download "$read2_fq" -o ${read2_root}.fq.gz
        echo "* Fastq file: '${read2_root}.fq.gz'"
    fi
    bam_root="${read1_root}_${read2_root}_bwa"

    bwa_ix_root=`dx describe "$bwa_index" --name`
    bwa_ix_root=${bwa_ix_root%.tar.gz}
    bwa_ix_root=${bwa_ix_root%.tgz}
    ref_id=${bwa_ix_root%_bwa_index}
    echo "* Downloading and extracting ${bwa_ix_root}.tgz file..."
    dx download "$bwa_index" -o ${bwa_ix_root}.tgz
    tar zxvf ${bwa_ix_root}.tgz
    echo "* Reference fasta: ${ref_id}.fa"

    echo "* Aligning with bwa..."
    set -x
    bwa aln -q 5 -l 32 -k 2 -t $nthreads ${ref_id} ${read1_root}.fq.gz > tmp_1.sai
    bwa aln -q 5 -l 32 -k 2 -t $nthreads ${ref_id} ${read2_root}.fq.gz > tmp_2.sai
    #bwa sampe ${ref_id} tmp_1.sai tmp_2.sai ${read1_root}.fq.gz ${read2_root}.fq.gz | samtools view -Shu - | samtools sort -m 5000000000 - ${prefix}
    bwa sampe ${ref_id} tmp_1.sai tmp_2.sai ${read1_root}.fq.gz ${read2_root}.fq.gz | samtools view -S -b /dev/stdin > tmp.bam
    set +x
    
    echo "* Sort bam..."
    set -x
    samtools sort -m 50000000000 tmp.bam ${bam_root}
    set +x

    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${bam_root}.bam > ${bam_root}_bam_qc.txt
    set +x
    #rm tmp.sai tmp.bam

    echo "* Prepare metadata json..."
    meta=`echo \"samtools_flagstats\": { `
    # 2142 + 0 in total (QC-passed reads + QC-failed reads)
    var=`grep "in total" ${bam_root}_bam_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta $var`
    # 0 + 0      duplicates
    var=`grep -w duplicates ${bam_root}_bam_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 204621 + 0 mapped (95.48%:-nan%)
    var=`grep -w "mapped" ${bam_root}_bam_qc.txt | head -1 | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep -w "mapped" ${bam_root}_bam_qc.txt | head -1 | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    meta=`echo $meta, \"mapped_pct\": \"$var\"`
    # 2142 + 0 paired in sequencing
    var=`grep "paired in sequencing" ${bam_root}_bam_qc.txt | awk '{printf "\"paired\": %d, \"paired_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 107149 + 0 read1
    var=`grep -w read1 ${bam_root}_bam_qc.txt | awk '{printf "\"read1\": %d, \"read1_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 107149 + 0 read2
    var=`grep -w read2 ${bam_root}_bam_qc.txt | awk '{printf "\"read2\": %d, \"read2_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 2046 + 0 properly paired (95.48%:-nan%)
    var=`grep "properly paired" ${bam_root}_bam_qc.txt | awk '{printf "\"paired\": %d, \"paired_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep "properly paired" ${bam_root}_bam_qc.txt | awk '{print $6}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    meta=`echo $meta, \"paired_pct\": \"$var\"`
    # 2046212 + 0 with itself and mate mapped
    # 0 + 0      singletons (0.00%:-nan%)
    var=`grep -w singletons ${bam_root}_bam_qc.txt | awk '{printf "\"singletons\": %d, \"singletons_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 0 + 0    with mate mapped to a different chr
    var=`grep "with mate mapped to a different chr" ${bam_root}_bam_qc.txt | grep -v "(mapQ>=5)" | awk '{printf "\"diff_chroms\": %d, \"diff_chroms_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 0 + 0    with mate mapped to a different chr (mapQ>=5)
    var=`grep "with mate mapped to a different chr (mapQ>=5)" ${bam_root}_bam_qc.txt | awk '{printf "\"diff_chroms_gt5\": %d, \"diff_chroms_gt5_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    meta=`echo $meta }`

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    details=`echo { $meta }`
    bam_bwa=$(dx upload ${bam_root}.bam --details "$details" --property QC="$meta" --brief)
    bam_bwa_qc=$(dx upload ${bam_root}_bam_qc.txt --brief)

    dx-jobutil-add-output bam_bwa "$bam_bwa" --class=file
    dx-jobutil-add-output bam_bwa_qc "$bam_bwa_qc" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
