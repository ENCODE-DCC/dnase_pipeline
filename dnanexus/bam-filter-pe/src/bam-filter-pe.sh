#!/bin/bash
# bam-filter-pe 0.0.1

main() {
    echo "* Installing phantompeakqualtools and dependencies (caTools, snow)..." 2>&1 | tee -a install.log
    set -x
    # phantompeakqualtools  : also resquires boost C libraries (on aws), boost C libraries (on aws) samtools (in resources/usr/bin)
    wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz -O phantompeakqualtools.tgz 2>&1 | tee -a install.log
    mkdir phantompeakqualtools
    tar -xzf phantompeakqualtools.tgz -C phantompeakqualtools --strip-components=1
    cd phantompeakqualtools
    echo "install.packages('bitops',dependencies=TRUE,repos='http://cran.cnr.berkeley.edu/')" > installPkgs.R
    echo "install.packages('http://cran.r-project.org/src/contrib/caTools_1.17.1.tar.gz',dependencies=TRUE)" > installPkgs.R
    echo "install.packages('http://cran.r-project.org/src/contrib/snow_0.3-13.tar.gz',dependencies=TRUE)" >> installPkgs.R
    echo "install.packages('spp_1.10.1.tar.gz',dependencies=TRUE)" >> installPkgs.R
    echo "* === installPkgs.R ===" 2>&1 | tee -a install.log
    cat installPkgs.R 2>&1 | tee -a install.log
    echo "* =====================" 2>&1 | tee -a install.log
    sudo Rscript installPkgs.R 2>&1 | tee -a install.log
    cd ..
    # additional executables in resources/usr/bin
    set +x
      
    echo "*****"
    echo "* Running: bam-filter-pe.sh v0.0.1"
    echo "* samtools: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* edwBamFilter: "`edwBamFilter 2>&1 | grep "edwBamFilter v" | awk '{print $2}'`
    #echo "*****"  The following could have been separate 'bam-stats' app
    echo "* edwBamStats: "`edwBamStats 2>&1 | grep "edwBamStats v" | awk '{print $2}'`
    echo "* phantompeakqualtools: "`grep Version phantompeakqualtools/README.txt | awk '{print $2}'`
    echo "* bedtools: "`bedtools --version 2>&1 | awk '{print $2}'`
    echo "*****"

    echo "* Value of bam_bwa: '$bam_bwa'"
    echo "* Value of map_thresh: '$map_thresh'"

    echo "* Download files..."
    bam_bwa_root=`dx describe "$bam_bwa" --name`
    bam_bwa_root=${bam_bwa_root%.bam}
    dx download "$bam_bwa" -o ${bam_bwa_root}.bam
    echo "* bam_bwa file: '${bam_bwa_root}.bam'"
    bam_filtered_root="${bam_bwa_root}_filtered"
    bam_no_chrM_root="${bam_bwa_root}_no_chrM"
    bam_sample_root="${bam_no_chrM_root}_${sample_size}_sample"

    # TODO: Need to make bam.bai?
    #echo "* Indexing bam..."
    #set -x
    #samtools index ${bam_bwa_root}.bam
    #set +x
    
    echo "* Filter on threashold..."
    set -x
    samtools view -F 1804 -f 2 -q ${map_thresh} -u ${bam_bwa_root}.bam | samtools sort -n - ${bam_filtered_root}_tmp
    samtools fixmate -O bam -r ${bam_filtered_root}_tmp.bam - | samtools view -F 1804 -f 2 -u - | samtools sort -m 50000000000  - ${bam_filtered_root}
    samtools index ${bam_filtered_root}.bam
    set +x

    echo "* Collect bam stats..."
    set -x
    samtools flagstat ${bam_filtered_root}.bam > ${bam_filtered_root}_qc.txt
    # Need samtools 20 !
    samtools stats ${bam_filtered_root}.bam > ${bam_filtered_root}_qc_full.txt
    cp ${bam_filtered_root}_qc.txt ${bam_filtered_root}_qc_full.txt
    set +x

    echo "* Prepare metadata for filtered bam..."
    meta=`echo \"samtools_flagstats\": { `
    # 2142 + 0 in total (QC-passed reads + QC-failed reads)
    var=`grep "in total" ${bam_filtered_root}_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta $var`
    # 0 + 0      duplicates
    var=`grep -w duplicates ${bam_filtered_root}_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 204621 + 0 mapped (95.48%:-nan%)
    var=`grep -w mapped ${bam_filtered_root}_qc.txt | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep -w mapped ${bam_filtered_root}_qc.txt | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    meta=`echo $meta, \"mapped_pct\": \"$var\"`
    # 2142 + 0 paired in sequencing
    var=`grep "paired in sequencing" ${bam_filtered_root}_qc.txt | awk '{printf "\"paired\": %d, \"paired_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 107149 + 0 read1
    var=`grep -w read1 ${bam_filtered_root}_qc.txt | awk '{printf "\"read1\": %d, \"read1_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 107149 + 0 read2
    var=`grep -w read2 ${bam_filtered_root}_qc.txt | awk '{printf "\"read2\": %d, \"read2_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 2046 + 0 properly paired (95.48%:-nan%)
    var=`grep "properly paired" ${bam_filtered_root}_qc.txt | awk '{printf "\"paired\": %d, \"paired_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep "properly paired" ${bam_filtered_root}_qc.txt | awk '{print $6}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    meta=`echo $meta, \"paired_pct\": \"$var\"`
    # 2046212 + 0 with itself and mate mapped
    # 0 + 0      singletons (0.00%:-nan%)
    var=`grep -w singletons ${bam_filtered_root}_qc.txt | awk '{printf "\"singletons\": %d, \"singletons_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 0 + 0    with mate mapped to a different chr
    var=`grep "with mate mapped to a different chr" ${bam_filtered_root}_qc.txt | grep -v "(mapQ>=5)" | awk '{printf "\"diff_chroms\": %d, \"diff_chroms_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 0 + 0    with mate mapped to a different chr (mapQ>=5)
    var=`grep "with mate mapped to a different chr (mapQ>=5)" ${bam_filtered_root}_qc.txt | grep -v "(mapQ>=5)" | awk '{printf "\"diff_chroms_gt5\": %d, \"diff_chroms_gt5_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    meta=`echo $meta }`
    qc_filtered=$meta

    # TODO: Collect QC from: ${bam_filtered_root}_qc_full.txt

    echo "* Filter out chrM..."
    set -x
    # Note the sort by name whih *may* be needed for the sampling/stats
    edwBamFilter -sponge -chrom=chrM ${bam_filtered_root}.bam ${bam_no_chrM_root}.bam  ## qc based on bam without chrm
    samtools sort-m 50000000000 -n ${bam_no_chrM_root}.bam ${bam_no_chrM_root}_byname ## for pbc usage
    samtools index ${bam_no_chrM_root}_byname.bam
    set +x

    #echo "*****"  The following could have been separate 'bam-stats' app

    echo "* Generating stats on $sample_size reads..."
    set -x
    edwBamStats -sampleBamSize=${sample_size} -u4mSize=${sample_size} -sampleBam=${bam_sample_root}.bam ${bam_no_chrM_root}_byname.bam bam_no_chrM_stats.txt
    samtools index ${bam_sample_root}.bam
    set +x

    echo "* Running spp..."
    set -x
    Rscript phantompeakqualtools/run_spp.R -x=-500:-1 -s=-500:5:1500 -rf -c=${bam_sample_root}.bam -out=${bam_sample_root}_spp.txt
    set +x

    echo "* Running pbc..."
    set -x
    #dac_pbc ${endedness} ${bam_sample_root}.bam > ${bam_sample_root}_pbc.txt
    ## m0 distinct reads; m1 location with exactly one read;
    ## m2 location with exactly 2 read; mt total reads
    # bamToBed -bedpe -i stdin \
    bamToBed -bedpe -i ${bam_sample_root}.bam | awk 'BEGIN{OFS="\t"} {print $1,$2,$4,$6,$9,$10}' | sort | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
        > ${bam_sample_root}_pbc.txt
    set +x

    # TODO: Collect QC from: ${bam_sample_root}_stats.txt, ${bam_sample_root}_spp.txt, ${bam_sample_root}_pbc.txt



    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    details=`echo { $qa_filtered }`
    bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "$details" --property QC="$qc_filtered" --brief)
    bam_no_chrM=$(dx upload ${bam_no_chrM_root}.bam --brief)
    bam_sample=$(dx upload ${bam_sample_root}.bam --brief)
    bam_filtered_qc=$(dx upload ${bam_filtered_root}_qc.txt --brief)
    bam_filtered_qc_full=$(dx upload ${bam_filtered_root}_qc_full.txt --brief)
    #echo "*****"  The following could have been separate 'bam-stats' app
    bam_sample_stats=$(dx upload ${bam_sample_root}_stats.txt --brief)
    bam_sample_spp=$(dx upload ${bam_sample_root}_spp.txt --brief)
    bam_sample_pbc=$(dx upload ${bam_sample_root}_pbc.txt --brief)

    dx-jobutil-add-output bam_filtered "$bam_filtered" --class=file
    dx-jobutil-add-output bam_no_chrM "$bam_no_chrM" --class=file
    dx-jobutil-add-output bam_filtered_qc "$bam_filtered_qc" --class=file
    dx-jobutil-add-output bam_filtered_qc_full "$bam_filtered_qc_full" --class=file
    #echo "*****"  The following could have been separate 'bam-stats' app
    dx-jobutil-add-output bam_sample_stats "$bam_sample_stats" --class=file
    dx-jobutil-add-output bam_sample_spp "$bam_sample_spp" --class=file
    dx-jobutil-add-output bam_sample_pbc "$bam_sample_pbc" --class=file

    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
