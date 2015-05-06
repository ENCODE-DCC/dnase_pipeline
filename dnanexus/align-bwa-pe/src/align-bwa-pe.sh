#!/bin/bash
# align-bwa-pe.sh

script_name="align-bwa-pe.sh"
script_ver="0.2.0"

main() {
    # Executable in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of read1_fq: '$read1_fq'"
    echo "* Value of read2_fq: '$read2_fq'"
    echo "* Value of bwa_index: '$bwa_index'"
    echo "* Value of nthreads: '$nthreads'"

    #echo "* Download files..."
    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!read1_fq[@]}
    do
        filename=`dx describe "${read1_fq[$ix]}" --name | cut -d'.' -f1`
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
        echo "* Downloading ${file_root}.fq.gz file..."
        dx download "${read1_fq[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Read1 fastq${concat} file: '${outfile_name}.fq.gz'"
    read1_root=${outfile_name}

    outfile_name=""
    concat=""
    rm -f concat.fq
    for ix in ${!read2_fq[@]}
    do
        filename=`dx describe "${read2_fq[$ix]}" --name | cut -d'.' -f1`
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
        echo "* Downloading ${file_root}.fq.gz file..."
        dx download "${read2_fq[$ix]}" -o - | gunzip >> concat.fq
    done
    mv concat.fq ${outfile_name}.fq
    echo "* Gzipping file..."
    gzip ${outfile_name}.fq
    echo "* Read2 fastq${concat} file: '${outfile_name}.fq.gz'"
    read2_root=${outfile_name}
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
    bwa sampe ${ref_id} tmp_1.sai tmp_2.sai ${read1_root}.fq.gz ${read2_root}.fq.gz \
                    | samtools view -Shu - | samtools sort -@ $nthreads -m 5G -f - tmp.sam
    samtools view -hb tmp.sam > ${bam_root}.bam
    samtools index ${bam_root}.bam
    #samtools view -H ${bam_root}.bam
    set +x
    
    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${bam_root}.bam > ${bam_root}_bam_qc.txt
    set +x
    #rm tmp.sai tmp.bam

    echo "* Prepare metadata json..."
    meta=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        meta=`qc_metrics.py -n samtools_flagstats -f ${bam_root}_bam_qc.txt`
    fi

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_bwa=$(dx upload ${bam_root}.bam --details "{ $meta }" --property QC="{ $meta }" --property SW="$versions" --brief)
    bam_bwa_qc=$(dx upload ${bam_root}_bam_qc.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_bwa "$bam_bwa" --class=file
    dx-jobutil-add-output bam_bwa_qc "$bam_bwa_qc" --class=file
    dx-jobutil-add-output metadata "{ $meta }" --class=string

    echo "* Finished."
}
