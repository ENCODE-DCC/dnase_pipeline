#!/bin/bash
# align-bwa-se 0.1.0

main() {
    # Executable in resources/usr/bin
    
    echo "*****"
    echo "* Running: align-bwa-se.sh v0.1.0"
    echo "* bwa version: "`bwa 2>&1 | grep Version | awk '{print $2}'`
    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "* Value of reads_fq: '$reads_fq'"
    echo "* Value of bwa_index: '$bwa_index'"
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
    if [ "${#reads_fq}" -gt "1" ]; then
        outfile_name="concat"
        rm -f concat.fq
        for ix in ${!reads_fq[@]}
        do
            filename=`dx describe "${reads_fq[$ix]}" --name | cut -d'.' -f1`
            file_root=${filename%.fastq.gz}
            file_root=${filename%.fq.gz}
            outfile_name="${file_root}_${outfile_name}"
            echo "* Downloading and concatenating ${file_root}.fq.gz file..."
            dx download "${reads_fq[$ix]}" -o - | gunzip >> concat.fq
        done
        mv concat.fq ${outfile_name}.fq
        echo "* Gzipping file..."
        gzip ${outfile_name}.fq
        echo "* Fastqs concatenated as: '${outfile_name}.fq.gz'"
        reads_root=${outfile_name}
    else
        reads_root=`dx describe "$reads_fq" --name`
        reads_root=${reads_root%.fastq.gz}
        reads_root=${reads_root%.fq.gz}
        echo "* Downloading ${reads_fq}.fq.gz file..."
        dx download "$reads_fq" -o ${reads_root}.fq.gz
        echo "* Fastq file: '${reads_root}.fq.gz'"
    fi
    bam_root="${reads_root}_bwa"

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
    bwa aln -q 5 -l 32 -k 2 -t $nthreads ${ref_id} ${reads_root}.fq.gz > tmp.sai
    #bwa samse ${ref_id} tmp.sai ${reads_root}.fq.gz | samtools view -Shu - | samtools sort -m 5000000000 - ${reads_root}
    bwa samse ${ref_id} tmp.sai ${reads_root}.fq.gz | samtools view -S -b /dev/stdin > tmp.bam
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
    # 5008100 + 0 in total (QC-passed reads + QC-failed reads)
    var=`grep "QC-passed reads" ${bam_root}_bam_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta $var`
    # 0 + 0 duplicates
    var=`grep -w duplicates ${bam_root}_bam_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 3421023 + 0 mapped (68.31%:-nan%)
    var=`grep -w "mapped" ${bam_root}_bam_qc.txt | head -1 | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep -w "mapped" ${bam_root}_bam_qc.txt | head -1 | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    meta=`echo $meta, \"mapped_pct\": \"$var\"`
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
