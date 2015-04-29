#!/bin/bash
# bam-filter-se.sh 0.1.0

script_name="bam-filter-se.sh"
script_ver="0.1.0"

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
    
    echo "*****"
    echo "* Running: $script_name: $script_ver"; versions=`echo "\"sw_versions\": { \"$script_name\": \"$script_ver\""`
    var=`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* samtools version: $var"; versions=`echo "$versions, \"samtools\": \"$var\""`
    var=`edwBamFilter 2>&1 | grep "edwBamFilter v" | awk '{print $2}'`
    echo "* edwBamFilter version: $var"; versions=`echo "$versions, \"edwBamFilter\": \"$var\""`
    var=`edwBamStats 2>&1 | grep "edwBamStats v" | awk '{print $2}'`
    echo "* edwBamStats version: $var"; versions=`echo "$versions, \"edwBamStats\": \"$var\""`
    #var=`R --version | grep "R version" | awk '{print $3,$4}'`
    #echo "* R version: $var"; versions=`echo "$versions, \"R\": \"$var\""`
    var=`Rscript --version 2>&1 | awk '{print $5,$6}'`
    echo "* Rscript version: $var"; versions=`echo "$versions, \"Rscript\": \"$var\""`
    var=`grep Version phantompeakqualtools/README.txt | awk '{print $2}'`
    echo "* phantompeakqualtools version: $var"; versions=`echo "$versions, \"phantompeakqualtools\": \"$var\""`
    var=`grep caTools_ phantompeakqualtools/install.log | head -1 | tr \_ " " | awk '{print $4}'`; var=`echo "${var%.tar*}"
    echo "* caTools version: $var"; versions=`echo "$versions, \"caTools\": \"$var\""`
    var=`grep snow_ phantompeakqualtools/install.log | head -1 | tr \_ " " | awk '{print $4}'`; var=`echo "${var%.tar*}"
    echo "* snow version: $var"; versions=`echo "$versions, \"snow\": \"$var\""`
    var=`grep spp_ phantompeakqualtools/installPkgs.R | tr \_ " " | awk '{print $2}'`; var=`echo "*${var%.tar*}"
    echo "* spp version: $var"; versions=`echo "$versions, \"spp\": \"$var\""`
    var=`gawk --version | grep Awk | awk '{print $3}'`
    echo "* gawk version: $var"; versions=`echo "$versions, \"gawk\": \"$var\""`
    var=`bedtools --version 2>&1 | awk '{print $2}'`
    echo "* bedtools version: $var"; versions=`echo "$versions, \"bedtools\": \"$var\""``
    versions=`echo $versions }`
    echo "*****"

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
    samtools view -b ${bam_filtered_root}.sam > ${bam_filtered_root}.bam
    samtools index ${bam_filtered_root}.bam
    set +x

    echo "* Collect filtered bam stats..."
    set -x
    samtools flagstat ${bam_filtered_root}.bam | tee ${bam_filtered_root}_qc.txt
    # cat ${bam_filtered_root}_qc.txt ## tee already did this
    samtools stats ${bam_filtered_root}.bam > ${bam_filtered_root}_qc_full.txt
    head -3 ${bam_filtered_root}_qc_full.txt
    grep ^SN ${bam_filtered_root}_qc_full.txt | cut -f 2- > ${bam_filtered_root}_qc_summary.txt
    cat ${bam_filtered_root}_qc_summary.txt
    set +x

    echo "* Prepare metadata for filtered bam..."
    meta=`echo \"samtools_flagstats\": { `
    # 2826233 + 0 in total (QC-passed reads + QC-failed reads)
    var=`grep "QC-passed reads" ${bam_filtered_root}_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta $var`
    # 0 + 0 duplicates
    var=`grep -w duplicates ${bam_filtered_root}_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    # 2826233 + 0 mapped (100.00%:-nan%)
    var=`grep -w "mapped" ${bam_filtered_root}_qc.txt | head -1 | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    meta=`echo $meta, $var`
    var=`grep -w "mapped" ${bam_filtered_root}_qc.txt | head -1 | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    meta=`echo $meta, \"mapped_pct\": \"$var\"`
    meta=`echo $meta }`
    qc_filtered=$meta

    meta=`echo \"samtools_stats\": { `
    # raw total sequences:	2826233
    var=`grep "^raw total sequences\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"raw total sequences\": %d", $4}'`
    meta=`echo $meta $var`
    # filtered sequences:	0
    var=`grep "^filtered sequences\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"filtered sequences\": %d", $3}'`
    meta=`echo $meta, $var`
    # sequences:	2826233
    var=`grep "^sequences\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"sequences\": %d", $2}'`
    meta=`echo $meta, $var`
    # is sorted:	1
    var=`grep "^is sorted\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"is sorted\": %d", $3}'`
    meta=`echo $meta, $var`
    # 1st fragments:	2826233
    var=`grep "^1st fragments\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"1st fragments\": %d", $3}'`
    meta=`echo $meta, $var`
    # last fragments:	0
    var=`grep "^last fragments\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"last fragments\": %d", $3}'`
    meta=`echo $meta, $var`
    # reads mapped:	2826233
    var=`grep "^reads mapped\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads mapped\": %d", $3}'`
    meta=`echo $meta, $var`
    # reads mapped and paired:	0	# paired-end technology bit set + both mates mapped
    var=`grep "^reads mapped and paired\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads mapped and paired\": %d", $5}'`
    meta=`echo $meta, $var`
    # reads unmapped:	0
    var=`grep "^reads unmapped\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads unmapped\": %d", $3}'`
    meta=`echo $meta, $var`
    # reads properly paired:	0	# proper-pair bit set
    var=`grep "^reads properly paired\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads properly paired\": %d", $4}'`
    meta=`echo $meta, $var`
    # reads paired:	0	# paired-end technology bit set
    var=`grep "^reads paired\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads paired\": %d", $3}'`
    meta=`echo $meta, $var`
    # reads duplicated:	0	# PCR or optical duplicate bit set
    var=`grep "^reads duplicated\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads duplicated\": %d", $3}'`
    meta=`echo $meta, $var`
    # reads MQ0:	0	$ mapped and MQ=0
    var=`grep "^reads MQ0\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads MQ0\": %d", $3}'`
    meta=`echo $meta, $var`
    # reads QC failed:	0
    var=`grep "^reads QC failed\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"reads QC failed\": %d", $4}'`
    meta=`echo $meta, $var`
    # non-primary alignments:	0
    var=`grep "^non-primary alignments\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"non-primary alignments\": %d", $3}'`
    meta=`echo $meta, $var`
    # total length:	101744388	# ignores clipping
    var=`grep "^total length\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"total length\": %d", $3}'`
    meta=`echo $meta, $var`
    # bases mapped:	101744388	# ignores clipping
    var=`grep "^bases mapped\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"bases mapped\": %d", $3}'`
    meta=`echo $meta, $var`
    # bases mapped (cigar):	101540364	# more accurate
    var=`grep "^bases mapped" ${bam_filtered_root}_qc_summary.txt | grep "cigar" | awk '{printf "\"bases mapped (cigar)\": %d", $4}'`
    meta=`echo $meta, $var`
    # bases trimmed:	0
    var=`grep "^bases trimmed\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"bases trimmed\": %d", $3}'`
    meta=`echo $meta, $var`
    # bases duplicated:	0
    var=`grep "^bases duplicated\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"bases duplicated\": %d", $3}'`
    meta=`echo $meta, $var`
    # mismatches:	474153	# from NM fields
    var=`grep "^mismatches\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"mismatches\": %d", $2}'`
    meta=`echo $meta, $var`
    # error rate:	4.669601e-03	# mismatches / bases mapped (cigar)
    var=`grep "^error rate\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"error rate\": %s", $3}'`
    meta=`echo $meta, $var`
    # average length:	36
    var=`grep "^average length\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"average length\": %s", $3}'`
    meta=`echo $meta, $var`
    # maximum length:	36
    var=`grep "^maximum length\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"maximum length\": %d", $3}'`
    meta=`echo $meta, $var`
    # average quality:	36.5
    var=`grep "^average quality\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"average quality\": %s", $3}'`
    meta=`echo $meta, $var`
    read_len=$var
    # insert size average:	0.0
    var=`grep "^insert size average\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"insert size average\": %s", $4}'`
    meta=`echo $meta, $var`
    # insert size standard deviation:	0.0
    var=`grep "^insert size standard deviation\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"insert size standard deviation\": %s", $5}'`
    meta=`echo $meta, $var`
    # inward oriented pairs:	0
    var=`grep "^inward oriented pairs\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"inward oriented pairs\": %d", $4}'`
    meta=`echo $meta, $var`
    # outward oriented pairs:	0
    var=`grep "^outward oriented pairs\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"outward oriented pairs\": %d", $4}'`
    meta=`echo $meta, $var`
    # pairs with other orientation:	0
    var=`grep "^pairs with other orientation\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"pairs with other orientation\": %d", $5}'`
    meta=`echo $meta, $var`
    # pairs on different chromosomes:	0
    var=`grep "^pairs on different chromosomes\:" ${bam_filtered_root}_qc_summary.txt | awk '{printf "\"pairs on different chromosomes\": %d", $5}'`
    meta=`echo $meta, $var`
    meta=`echo $meta }`
    qc_filtered=`echo $qc_filtered, $meta`

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
    cat ${bam_sample_root}_stats.txt
    samtools index ${bam_sample_root}.bam
    cat ${bam_sample_root}_stats.txt
    set +x

    echo "* Running spp..."
    set -x
    # awk didn't work, so use gawk and pre-create the tagAlign 
    samtools view -F 0x0204 -o - ${bam_sample_root}.bam | \
       gawk 'BEGIN{OFS="\t"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),"N","1000","-"} else {print $3,($4-1),($4-1+length($10)),"N","1000","+"} }' \
       | gzip -c > ${bam_sample_root}.tagAlign.gz
    Rscript phantompeakqualtools/run_spp.R -x=-500:-1 -s=-500:5:1500 -rf -c=${bam_sample_root}.tagAlign.gz \
                                        -out=${bam_sample_root}_spp.txt | tee ${bam_sample_root}_spp_out.txt
    #cat ${bam_sample_root}_spp.txt
    #ENCFF001RUR_rep1_1_bwa_no_chrM_15000000_sample.bam	2286836	40,55	0.086978740217396,0.0864436517524779	40	0.08697874	1500	0.03001068	2.898259	1	1
    touch ${bam_sample_root}_spp.txt
    set +x

    echo "* Running pbc..."
    set -x
    bedtools bamtobed -i ${bam_sample_root}.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | sort | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' \
        > ${bam_sample_root}_pbc.txt
    cat ${bam_sample_root}_pbc.txt 
    # Interpretation??  mt=read_count m0=locations_mapped m1=mapped_by_one m2=mapped_by_two m0/mt=locations_per_read m1/m0=ratio_1_over_locations m1/m2=ratio_1_over_2
    # 2286836	2219898	2176268	37175	0.970729	0.980346	58.541170
    set +x

    echo "* Prepare metadata for sampled bam..."
    meta=`echo \"edwBamStats\": { `
    # alignedBy bwa
    var=`grep "^alignedBy" ${bam_sample_root}_stats.txt | awk '{printf "\"alignedBy\": \"%s\"", $2}'`
    meta=`echo $meta $var`
    # isPaired 0
    var=`grep "^isPaired" ${bam_sample_root}_stats.txt | awk '{printf "\"isPaired\": %d", $2}'`
    meta=`echo $meta, $var`
    # isSortedByTarget 1
    var=`grep "^isSortedByTarget" ${bam_sample_root}_stats.txt | awk '{printf "\"isSortedByTarget\": %d", $2}'`
    meta=`echo $meta, $var`
    # readCount 2286836
    var=`grep "^readCount" ${bam_sample_root}_stats.txt | awk '{printf "\"readCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # readBaseCount 82326096
    var=`grep "^readBaseCount" ${bam_sample_root}_stats.txt | awk '{printf "\"readBaseCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # mappedCount 2286836
    var=`grep "^mappedCount" ${bam_sample_root}_stats.txt | awk '{printf "\"mappedCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # uniqueMappedCount 2285489
    var=`grep "^uniqueMappedCount" ${bam_sample_root}_stats.txt | awk '{printf "\"uniqueMappedCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeMean 36
    var=`grep "^readSizeMean" ${bam_sample_root}_stats.txt | awk '{printf "\"readSizeMean\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeStd 0
    var=`grep "^readSizeStd" ${bam_sample_root}_stats.txt | awk '{printf "\"readSizeStd\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeMin 36
    var=`grep "^readSizeMin" ${bam_sample_root}_stats.txt | awk '{printf "\"readSizeMin\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeMax 36
    var=`grep "^readSizeMax" ${bam_sample_root}_stats.txt | awk '{printf "\"readSizeMax\": %d", $2}'`
    meta=`echo $meta, $var`
    # u4mReadCount 2285489
    var=`grep "^u4mReadCount" ${bam_sample_root}_stats.txt | awk '{printf "\"u4mReadCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # u4mUniquePos 2213397
    var=`grep "^u4mUniquePos" ${bam_sample_root}_stats.txt | awk '{printf "\"u4mUniquePos\": %d", $2}'`
    meta=`echo $meta, $var`
    # u4mUniqueRatio 0.968457
    var=`grep "^u4mUniqueRatio" ${bam_sample_root}_stats.txt | awk '{printf "\"u4mUniqueRatio\": %s", $2}'`
    meta=`echo $meta, $var`
    # targetSeqCount 25
    var=`grep "^targetSeqCount" ${bam_sample_root}_stats.txt | awk '{printf "\"targetSeqCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # targetBaseCount 3095693983
    var=`grep "^targetBaseCount" ${bam_sample_root}_stats.txt | awk '{printf "\"targetBaseCount\": %d", $2}'`
    meta=`echo $meta, $var`
    meta=`echo $meta }`
    qc_sampled=$meta
    
    # Note that all values in ${bam_sample_root}_spp.txt are found in ${bam_sample_root}_spp_out.txt
    meta=`echo \"phantompeaktools_spp\": { `
    # strandshift(min): -500 
    var=`grep "^strandshift" ${bam_sample_root}_spp_out.txt | grep "min" | awk '{printf "\"strandshift(min)\": %s", $2}'`
    meta=`echo $meta $var`
    # strandshift(step): 5 
    var=`grep "^strandshift" ${bam_sample_root}_spp_out.txt | grep "step" | awk '{printf "\"strandshift(step)\": %s", $2}'`
    meta=`echo $meta, $var`
    # strandshift(max) 1500 
    var=`grep "^strandshift" ${bam_sample_root}_spp_out.txt | grep "max" | awk '{printf "\"strandshift(max)\": %s", $2}'`
    meta=`echo $meta, $var`
    # exclusion(min): -500 
    var=`grep "^exclusion" ${bam_sample_root}_spp_out.txt | grep "min" | awk '{printf "\"exclusion(min)\": %s", $2}'`
    meta=`echo $meta, $var`
    # exclusion(max): -1 
    var=`grep "^exclusion" ${bam_sample_root}_spp_out.txt | grep "max" | awk '{printf "\"exclusion(max)\": %s", $2}'`
    meta=`echo $meta, $var`
    # FDR threshold: 0.01 
    var=`grep "^FDR threshold\:" ${bam_sample_root}_spp_out.txt | awk '{printf "\"FDR threshold\": %s", $3}'`
    meta=`echo $meta, $var`
    # done. read 2286836 fragments
    var=`grep "^done\. read" ${bam_sample_root}_spp_out.txt | awk '{printf "\"read fragments\": %s", $3}'`
    meta=`echo $meta, $var`
    # ChIP data read length 36 
    var=`grep "^ChIP data read length" ${bam_sample_root}_spp_out.txt | awk '{printf "\"ChIP data read length\": %s", $5}'`
    meta=`echo $meta, $var`
    # Minimum cross-correlation value 0.03001068 
    var=`grep "^Minimum cross-correlation value" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Minimum cross-correlation value\": %s", $4}'`
    meta=`echo $meta, $var`
    # Minimum cross-correlation shift 1500 
    var=`grep "^Minimum cross-correlation shift" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Minimum cross-correlation shift\": %s", $4}'`
    meta=`echo $meta, $var`
    # Top 3 cross-correlation values 0.086978740217396,0.0864436517524779 
    var=`grep "^Top 3 cross-correlation values" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Top 3 cross-correlation values\": [ %s ]", $5}'`
    meta=`echo $meta, $var`
    # Top 3 estimates for fragment length 40,55 
    var=`grep "^Top 3 estimates for fragment length" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Top 3 estimates for fragment length\": [ %s ]", $7}'`
    meta=`echo $meta, $var`
    # Window half size 265 
    var=`grep "^Window half size" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Window half size\": %s", $4}'`
    meta=`echo $meta, $var`
    # Phantom peak location 40 
    var=`grep "^Phantom peak location" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Phantom peak location\": %s", $4}'`
    meta=`echo $meta, $var`
    # Phantom peak Correlation 0.08697874 
    var=`grep "^Phantom peak Correlation" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Phantom peak Correlation\": %s", $4}'`
    meta=`echo $meta, $var`
    # Normalized Strand cross-correlation coefficient (NSC) 2.898259 
    var=`grep "^Normalized Strand cross-correlation" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Normalized Strand cross-correlation coefficient (NSC)\": %s", $6}'`
    meta=`echo $meta, $var`
    # Relative Strand cross-correlation Coefficient (RSC) 1 
    var=`grep "^Relative Strand cross-correlation" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Relative Strand cross-correlation Coefficient (RSC)\": %s", $6}'`
    meta=`echo $meta, $var`
    # Phantom Peak Quality Tag 1 
    var=`grep "^Phantom Peak Quality Tag" ${bam_sample_root}_spp_out.txt | awk '{printf "\"Phantom Peak Quality Tag\": %s", $5}'`
    meta=`echo $meta, $var`
    meta=`echo $meta }`
    qc_sampled=`echo $qc_sampled, $meta`

    #cat ${bam_sample_root}_pbc.txt 
    # Interpretation??  mt="reads" m0="locations" m1="mapped by 1 read" m2="mapped by 2 reads" m0/mt="locations per read" m1/m0="proportion of 1 read locations" m1/m2="ratio: 1 read over 2 read locations"
    # 2286836	2219898	2176268	37175	0.970729	0.980346	58.541170
    meta=`echo \"pbc\": { `
    var=`cat ${bam_sample_root}_pbc.txt | awk '{printf "\"reads\": %s, \"locations\": %s, \"mapped by 1 read\": %s, \"mapped by 2 reads\": %s, \"locations per read\": %s, \"proportion of 1 read locations\": %s, \"ratio: 1 read over 2 read locations\": %s", $1,$2,$3,$4,$5,$6,$7}'`
    meta=`echo $meta $var`
    meta=`echo $meta }`
    qc_sampled=`echo $qc_sampled, $meta`
    
    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "{ $qc_filtered }" --property QC="$qc_filtered" --property SW="$versions" --brief)
    bam_no_chrM=$(dx upload ${bam_no_chrM_root}.bam --details "{ $read_len }" --property SW="$versions" --brief)
    bam_sample=$(dx upload ${bam_sample_root}.bam --details "{ $qc_sampled }" --property QC="$qc_sampled" --property SW="$versions" --brief)
    bam_filtered_qc=$(dx upload ${bam_filtered_root}_qc.txt --property SW="$versions" --brief)
    bam_filtered_qc_full=$(dx upload ${bam_filtered_root}_qc_full.txt --property SW="$versions" --brief)
    bam_sample_stats=$(dx upload ${bam_sample_root}_stats.txt --property SW="$versions" --brief)
    bam_sample_spp=$(dx upload ${bam_sample_root}_spp.txt --property SW="$versions" --brief)
    bam_sample_pbc=$(dx upload ${bam_sample_root}_pbc.txt --property SW="$versions" --brief)

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
