#!/bin/bash
# dnase-eval-bam-pe.sh - Evaluates sample of (paired-end) bam for the ENCODE DNase-seq pipeline.

main() {
    echo "* Installing phantompeakqualtools, caTools, snow and spp..." 2>&1 | tee -a install.log
    set -x
    # phantompeakqualtools  : also resquires boost C libraries (on aws), boost C libraries (on aws) samtools (in resources/usr/bin)
    wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz -O phantompeakqualtools.tgz >> install.log 2>&1
    mkdir phantompeakqualtools
    tar -xzf phantompeakqualtools.tgz -C phantompeakqualtools --strip-components=1
    cd phantompeakqualtools
    # By not having caTools and snow in execDepends, we can at least show which versions are used
    echo "install.packages(c('caTools','snow'),dependencies=TRUE,repos='http://cran.cnr.berkeley.edu/')" > installPkgs.R
    echo "install.packages('spp_1.10.1.tar.gz')" >> installPkgs.R
    cat installPkgs.R
    sudo Rscript installPkgs.R >> install.log 2>&1
    grep DONE install.log
    set +x
    done_count=`grep -c DONE install.log`
    if [ "$done_count" != "6" ]; then
        set -x
        cat install.log
        set -x
    fi
    set -x
    cd ..
    # additional executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi
 
    echo "* Value of bam_filtered: '$bam_filtered'"
    echo "* Value of sample_size: '$sample_size'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    # expecting *_bwa_biorep_filtered.bam
    bam_input_root=`dx describe "$bam_filtered" --name`
    # Better to leave the whole suffix!
    bam_input_root=${bam_input_root%.bam}
    #bam_input_root=${bam_input_root%_sized}
    #bam_input_root=${bam_input_root%_filtered}
    #bam_input_root=${bam_input_root%_biorep}
    #bam_input_root=${bam_input_root%_techrep}
    #bam_input_root=${bam_input_root%_bwa}
    dx download "$bam_filtered" -o ${bam_input_root}.bam
    echo "* bam file: '${bam_input_root}.bam'"
    bam_no_chrM_root="${bam_input_root}_no_chrM"
    bam_sample_root="${bam_input_root}_${sample_size}_sample"
    
    echo "* Filter out chrM..."
    # Note the sort by name which is needed for proper pe sampling
    set -x
    edwBamFilter -sponge -chrom=chrM ${bam_input_root}.bam ${bam_no_chrM_root}.bam  ## qc based on bam without chrm
    samtools sort -@ $nthreads -m 6G -n -f ${bam_no_chrM_root}.bam ${bam_no_chrM_root}_byname.sam ## for pbc usage
    samtools view -hb ${bam_no_chrM_root}_byname.sam > ${bam_no_chrM_root}_byname.bam
    samtools index ${bam_no_chrM_root}_byname.bam
    rm *.sam
    set +x

    echo "* Generating stats on $sample_size reads..."
    set -x
    edwBamStats -sampleBamSize=${sample_size} -u4mSize=${sample_size} -sampleBam=${bam_sample_root}.bam \
                                              ${bam_no_chrM_root}_byname.bam ${bam_no_chrM_root}_sampling_edwBamStats.txt
    edwBamStats ${bam_sample_root}.bam ${bam_sample_root}_edwBamStats.txt
    samtools index ${bam_sample_root}.bam
    set +x

    echo "* Running spp..."
    set -x
    # awk didn't work, so use gawk and pre-create the tagAlign 
    samtools view -F 0x0204 -o - ${bam_sample_root}.bam | \
       gawk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' \
       | gzip -c > ${bam_sample_root}.tagAlign.gz
    Rscript phantompeakqualtools/run_spp.R -x=-500:-1 -s=-500:5:1500 -rf -c=${bam_sample_root}.tagAlign.gz \
                                        -out=${bam_sample_root}_spp.txt > ${bam_sample_root}_spp_out.txt
    touch ${bam_sample_root}_spp.txt
    set +x

    echo "* Running pbc..."
    set -x
    bedtools bamtobed -i ${bam_sample_root}.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
        > ${bam_sample_root}_pbc.txt
    set +x

    echo "* Prepare metadata for sampled bam..."
    qc_sampled=''
    reads_sampled=$sample_size
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_sampled=`qc_metrics.py -n edwBamStats -f ${bam_sample_root}_edwBamStats.txt`
        reads_sampled=`qc_metrics.py -n edwBamStats -f ${bam_sample_root}_edwBamStats.txt -k readCount`
        # Note that all values in ${bam_sample_root}_spp.txt are found in ${bam_sample_root}_spp_out.txt
        meta=`qc_metrics.py -n phantompeaktools_spp -f ${bam_sample_root}_spp_out.txt`
        qc_sampled=`echo $qc_sampled, $meta`
        meta=`qc_metrics.py -n pbc -f ${bam_sample_root}_pbc.txt`
        qc_sampled=`echo $qc_sampled, $meta`
    fi
    # All qc to one file per target file:
    echo "===== edwBamStats ====="           > ${bam_sample_root}_qc.txt
    cat ${bam_sample_root}_edwBamStats.txt  >> ${bam_sample_root}_qc.txt
    echo " "                                >> ${bam_sample_root}_qc.txt
    echo "===== phantompeaktools spp =====" >> ${bam_sample_root}_qc.txt
    cat ${bam_sample_root}_spp_out.txt      >> ${bam_sample_root}_qc.txt
    echo " "                                >> ${bam_sample_root}_qc.txt
    echo "===== bedtools pbc  ====="        >> ${bam_sample_root}_qc.txt
    cat ${bam_sample_root}_pbc.txt          >> ${bam_sample_root}_qc.txt
        
    echo "* Upload results..."
    bam_sample=$(dx upload ${bam_sample_root}.bam --details "{ $qc_sampled }" --property SW="$versions" \
                           --property sampled_reads="$reads_sampled" --property read_length="$read_len" --brief)
    bam_sample_qc=$(dx upload ${bam_sample_root}_qc.txt   --details "{ $qc_sampled }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_sample "$bam_sample" --class=file
    dx-jobutil-add-output bam_sample_qc "$bam_sample_qc" --class=file

    dx-jobutil-add-output sampled_reads "$reads_sampled" --class=string
    dx-jobutil-add-output metadata "{ $qc_sampled }" --class=string

    echo "* Finished."
}
