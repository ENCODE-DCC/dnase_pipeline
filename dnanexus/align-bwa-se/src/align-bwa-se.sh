#!/bin/bash
# align-bwa-se 0.0.1

main() {
    # Executable in resources/usr/bin
    
    echo "*****"
    echo "* Running: align-bwa-se.sh v0.0.1"
    echo "* bwa: "`bwa 2>&1 | grep Version | awk '{print $2}'`
    echo "* samtools: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "* Value of reads_fq: '$reads_fq'"
    echo "* Value of bwa_index: '$bwa_index'"
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
    if [ ${#reads_fq} gt 1 ]; then
        outfile_name="_concat"
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
    bwa_ix_root=${bwa_ix_root%.fasta.gz}
    bwa_ix_root=${bwa_ix_root.fa.gz}
    echo "* Downloading and unzipping ${bwa_ix_root}.fa.gz file..."
    dx download "$bwa_index" -o bwa_index.fa.gz
    gunzip bwa_index.fa.gz
    echo "* bwa index file: 'bwa_index.fa'"

    echo "* Aligning with bwa..."
    set -x
    bwa aln -q 5 -l 32 -k 2 -t $nthreads bwa_index.fa ${reads_root}.fq.gz > tmp.sai
    #bwa samse bwa_index.fa tmp.sai ${reads_root}.fq.gz | samtools view -Shu - | samtools sort -m 5000000000 - ${reads_root}
    bwa samse bwa_index.fa.gz tmp.sai ${reads_root}.fq.gz | samtools view -S -b /dev/stdin > tmp.bam
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
    # 2142994 + 0 in total (QC-passed reads + QC-failed reads)
    var=`grep "QC-passed reads" ${bam_root}_bam_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta $var`
    # 0 + 0 duplicates
    var=`grep -w duplicates ${bam_root}_bam_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 2046212 + 0 mapped (95.48%:-nan%)
    var=`grep -w mapped ${bam_root}_bam_qc.txt | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep -w mapped ${bam_root}_bam_qc.txt | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
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
