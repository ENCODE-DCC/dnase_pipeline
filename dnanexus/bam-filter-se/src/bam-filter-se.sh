#!/bin/bash
# bam-filter-se.sh

script_name="bam-filter-se.sh"
script_ver="0.2.0"

main() {
    echo "* Installing phantompeakqualtools, caTools, snow and spp..." 2>&1 | tee -a install.log
    set -x
    # phantompeakqualtools  : also resquires boost C libraries (on aws), boost C libraries (on aws) samtools (in resources/usr/bin)
    ls -l
    wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz -O phantompeakqualtools.tgz >> install.log
    mkdir phantompeakqualtools
    tar -xzf phantompeakqualtools.tgz -C phantompeakqualtools --strip-components=1
    cd phantompeakqualtools
    # By not having caTools and snow in execDepends, we can at least show which versions are used
    echo "install.packages(c('caTools','snow'),dependencies=TRUE,repos='http://cran.cnr.berkeley.edu/')" > installPkgs.R
    echo "install.packages('spp_1.10.1.tar.gz')" >> installPkgs.R
    sudo Rscript installPkgs.R >> install.log 2>&1
    cd ..
    # additional executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --applet $script_name --appver $script_ver`
    fi

    echo "* Value of bam_bwa: '$bam_bwa'"
    echo "* Value of map_thresh: '$map_thresh'"
    echo "* Value of sample_size: '$sample_size'"
    echo "* Value of nthreads: '$nthreads'"

    echo "* Download files..."
    bam_bwa_root=`dx describe "$bam_bwa" --name`
    bam_bwa_root=${bam_bwa_root%.bam}
    dx download "$bam_bwa" -o ${bam_bwa_root}.bam
    echo "* bam_bwa file: '${bam_bwa_root}.bam'"
    bam_filtered_root="${bam_bwa_root}_filtered"
    bam_no_chrM_root="${bam_bwa_root}_no_chrM"
    bam_sample_root="${bam_no_chrM_root}_${sample_size}_sample"
    
    echo "* Filter on threashold..."
    set -x
    samtools view -F 1804 -q ${map_thresh} -u ${bam_bwa_root}.bam | \
            samtools sort -@ $nthreads -m 5G -f - ${bam_filtered_root}.sam
    samtools view -hb ${bam_filtered_root}.sam > ${bam_filtered_root}.bam
    samtools index ${bam_filtered_root}.bam
    set +x

    echo "* Collect filtered bam stats..."
    set -x
    samtools flagstat ${bam_filtered_root}.bam > ${bam_filtered_root}_qc.txt
    samtools stats ${bam_filtered_root}.bam > ${bam_filtered_root}_qc_full.txt
    head -3 ${bam_filtered_root}_qc_full.txt
    grep ^SN ${bam_filtered_root}_qc_full.txt | cut -f 2- > ${bam_filtered_root}_qc_summary.txt
    set +x

    echo "* Prepare metadata for filtered bam..."
    qc_filtered=''
    read_len=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_filtered=`qc_metrics.py -n samtools_flagstats -f ${bam_filtered_root}_qc.txt`
        meta=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_filtered_root}_qc_summary.txt`
        read_len=`qc_metrics.py -n samtools_stats -d ':' -f ${bam_filtered_root}_qc_summary.txt --keypair "average length"`
        qc_filtered=`echo $qc_filtered, $meta`
    fi

    echo "* Filter out chrM..."
    # Note, unlike in pe there is no sort by name
    set -x
    edwBamFilter -sponge -chrom=chrM ${bam_filtered_root}.bam ${bam_no_chrM_root}.bam  ## qc based on bam without chrm
    samtools index ${bam_no_chrM_root}.bam
    set +x

    #echo "*****"  The following could have been separate 'bam-stats' app

    echo "* Generating stats on $sample_size reads..."
    set -x
    edwBamStats -sampleBamSize=${sample_size} -u4mSize=${sample_size} -sampleBam=${bam_sample_root}.bam \
                                                        ${bam_no_chrM_root}.bam ${bam_sample_root}_stats.txt
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
    set +x

    echo "* Running pbc..."
    set -x
    bedtools bamtobed -i ${bam_sample_root}.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
        > ${bam_sample_root}_pbc.txt
    set +x

    echo "* Prepare metadata for sampled bam..."
    qc_sampled=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_sampled=`qc_metrics.py -n edwBamStats -f ${bam_sample_root}_stats.txt`
        # Note that all values in ${bam_sample_root}_spp.txt are found in ${bam_sample_root}_spp_out.txt
        meta=`qc_metrics.py -n phantompeaktools_spp -f ${bam_sample_root}_spp_out.txt`
        qc_sampled=`echo $qc_sampled, $meta, $read_len`
        meta=`qc_metrics.py -n pbc -f ${bam_sample_root}_pbc.txt`
        qc_sampled=`echo $qc_sampled, $meta`
    fi
    
    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "{ $qc_filtered }" --property QC="{ $qc_filtered }" --property SW="$versions" --brief)
    bam_no_chrM=$(dx upload ${bam_no_chrM_root}.bam   --details "{ $read_len }"                                     --property SW="$versions" --brief)
    bam_sample=$(dx upload ${bam_sample_root}.bam     --details "{ $qc_sampled }"  --property QC="{ $qc_sampled }"  --property SW="$versions" --brief)
    bam_filtered_qc=$(dx upload ${bam_filtered_root}_qc.txt           --property SW="$versions" --brief)
    bam_filtered_qc_full=$(dx upload ${bam_filtered_root}_qc_full.txt --property SW="$versions" --brief)
    bam_sample_stats=$(dx upload ${bam_sample_root}_stats.txt         --property SW="$versions" --brief)
    bam_sample_spp=$(dx upload ${bam_sample_root}_spp.txt             --property SW="$versions" --brief)
    #bam_sample_spp=$(dx upload ${bam_sample_root}_spp_out.txt         --property SW="$versions" --brief)
    bam_sample_pbc=$(dx upload ${bam_sample_root}_pbc.txt             --property SW="$versions" --brief)

    dx-jobutil-add-output bam_filtered "$bam_filtered" --class=file
    dx-jobutil-add-output bam_no_chrM "$bam_no_chrM" --class=file
    dx-jobutil-add-output bam_sample "$bam_sample" --class=file
    dx-jobutil-add-output bam_filtered_qc "$bam_filtered_qc" --class=file
    dx-jobutil-add-output bam_filtered_qc_full "$bam_filtered_qc_full" --class=file
    dx-jobutil-add-output bam_sample_stats "$bam_sample_stats" --class=file
    dx-jobutil-add-output bam_sample_spp "$bam_sample_spp" --class=file
    dx-jobutil-add-output bam_sample_pbc "$bam_sample_pbc" --class=file

    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
