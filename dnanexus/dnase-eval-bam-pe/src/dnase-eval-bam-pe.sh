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
    
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_eval_bam_pe.sh ${bam_input_root}.bam $sample_size $nthreads
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    bam_sample_root="${bam_input_root}_${sample_size}_sample"

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
