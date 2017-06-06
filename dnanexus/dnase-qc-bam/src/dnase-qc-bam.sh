#!/bin/bash
# dnase-eval-bam.sh - Evaluates sample of bam for the ENCODE DNase-seq pipeline.

main() {
    echo "* Installing phantompeakqualtools, caTools, snow and spp..." 2>&1 | tee -a install.log
    set -x
    # phantompeakqualtools  : also resquires boost C libraries (on aws), boost C libraries (on aws) samtools (in resources/usr/bin)
    wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz -O phantompeakqualtools.tgz | tee -a install.log 2>&1
    mkdir phantompeakqualtools
    #tar -xzf phantompeakqualtools.tgz -C phantompeakqualtools --strip-components=1
    tar -xzf phantompeakqualtools.tgz
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
    cd ..
    exe_dir="`pwd`"
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz -O gsl.tgz >> install.log 2>&1
    mkdir gsl
    tar -xzf gsl.tgz -C gsl --strip-components=1
    cd gsl
    ./configure >> install.log 2>&1
    make > install.log 2>&1
    sudo make install >> install.log 2>&1
    gsl-config --libs > install.log 2>&1
    export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
    cd ..
    #wget https://github.com/rthurman/hotspot/archive/v4.1.0.tar.gz -O hotspot.tgz >> install.log 2>&1
    wget https://github.com/StamLab/hotspot/archive/v4.1.1.tar.gz -O hotspot.tgz >> install.log 2>&1
    mkdir hotspot
    tar -xzf hotspot.tgz -C hotspot --strip-components=1
    cd hotspot/hotspot-distr/hotspot-deploy
    make >> install.log 2>&1
    # Can either put bin in path or copy contents of bin to /usr/bin
    export PATH=${exe_dir}/hotspot/hotspot-distr/hotspot-deploy/bin:${PATH}
    cd ../../../
    set +x; 
    # additional executables in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi
 
    echo "* Value of bam_filtered:     '$bam_filtered'"
    echo "* Value of pe_or_se:         '$pe_or_se'"
    echo "* Value of sample_size:      '$sample_size'"
    echo "* Value of hotspot_mappable: '$hotspot_mappable'"
    #echo "* Value of read_length:     '$read_length'"
    echo "* Value of genome:           '$genome'"
    echo "* Value of nthreads:         '$nthreads'"

    echo "* Download files..."
    # expecting *_bwa_biorep_filtered.bam
    bam_input_root=`dx describe "$bam_filtered" --name`
    # Better to leave the whole suffix!
    bam_input_root=${bam_input_root%.bam}
    dx download "$bam_filtered" -o ${bam_input_root}.bam
    echo "* bam file: '${bam_input_root}.bam'"
    
    #dx download "$chrom_sizes" -o chrom.sizes
    mappable_archive=`dx describe "$hotspot_mappable" --name`
    dx download "$hotspot_mappable"

    # check read_length
    read_length=36   # Only supporting read_length 36 at this time
    #if [ -f /usr/bin/parse_property.py ]; then 
    #    read_len=`parse_property.py -f "$bam_file" -p "read_length" --quiet`
    #    if [ "$read_len" != "" ] && [ "$read_length" -ne "$read_len" ]; then
    #        if [ "$read_len" == "32" ] || [ "$read_len" == "36" ] || [ "$read_len" == "40" ] || [ "$read_len" == "50" ] \
    #        || [ "$read_len" == "58" ] || [ "$read_len" == "72" ] || [ "$read_len" == "76" ] || [ "$read_len" == "100" ]; then
    #            echo "* NOTE: Read length ($read_length) does not match discovered read size ($read_len). Using $read_len."
    #            read_length=$read_len
    #        else
    #            echo "* WARNING: Read length ($read_length) does not match discovered read size ($read_len)."
    #        fi
    #    fi
    #fi

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_qc_bam.sh ${bam_input_root}.bam $sample_size $nthreads $pe_or_se $genome $mappable_archive $read_length hotspot/
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
        meta=`qc_metrics.py -n pbc_spp --string "{ $qc_sampled }"`
        qc_sampled=`echo $qc_sampled, $meta`
        meta=`qc_metrics.py -n hotspot1 -f ${bam_sample_root}_hotspot1_qc.txt`
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
    echo -e "Sampled Reads\tDistinct Locations Mapped\tSingle-read Locations\tMulti-read Locations\tNRF (Non-Redundant Fraction)=Distinct Locations/Sample Reads\tPBC1=Single-read Locations/Distinct Locations\tPBC2=Single-read Locations/Multi-read Locations" >> ${bam_sample_root}_qc.txt
    cat ${bam_sample_root}_pbc.txt          >> ${bam_sample_root}_qc.txt
    echo "===== HotSpot1 ====="             >> ${bam_sample_root}_qc.txt
    cat ${bam_sample_root}_hotspot1_qc.txt  >> ${bam_sample_root}_qc.txt
        
    echo "* Upload results..."
    bam_sample=$(dx upload ${bam_sample_root}.bam --details "{ $qc_sampled }" --property SW="$versions" \
                            --property pe_or_se=$pe_or_se --property sampled_reads="$reads_sampled" \
                            --property read_length="$read_len" --brief)
    bam_sample_qc=$(dx upload ${bam_sample_root}_qc.txt   --details "{ $qc_sampled }" --property SW="$versions" --brief)

    dx-jobutil-add-output bam_sample "$bam_sample" --class=file
    dx-jobutil-add-output bam_sample_qc "$bam_sample_qc" --class=file

    dx-jobutil-add-output sampled_reads "$reads_sampled" --class=string
    dx-jobutil-add-output metadata "{ $qc_sampled }" --class=string

    echo "* Finished."
}
