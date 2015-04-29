#!/bin/bash
# sample-hotspots

script_name="sample-hotspots.sh"
script_ver="0.1.0"

main() {
    echo "* Installing hotspot and dependencies (gsl)..." 2>&1 | tee -a install.log
    exe_dir="`pwd`"
    set -x
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
    wget https://github.com/rthurman/hotspot/archive/v4.1.0.tar.gz -O hotspot.tgz >> install.log 2>&1
    mkdir hotspot
    tar -xzf hotspot.tgz -C hotspot --strip-components=1
    cd hotspot/hotspot-distr/hotspot-deploy
    make >> install.log 2>&1
    # Can either put bin in path or copy contents of bin to /usr/bin
    export PATH=${exe_dir}/hotspot/hotspot-distr/hotspot-deploy/bin:${PATH}
    cd ../../../
    set +x; 
    # additional executables in resources/usr/bin
    
    echo "*****"
    echo "* Running: $script_name: $script_ver"; versions=`echo "\"sw_versions\": { \"$script_name\": \"$script_ver\""`
    var=`edwBamStats 2>&1 | grep "edwBamStats v" | awk '{print $2}'`
    echo "* edwBamStats version: $var"; versions=`echo "$versions, \"edwBamStats\": \"$var\""`
    var=`hotspot 2>&1 | grep HotSpot | awk '{print $1}'`
    echo "* hotspot version: $var"; versions=`echo "$versions, \"hotspot\": \"$var\""`
    var=`python2.7 /usr/bin/hotspot.py -h | grep Version | awk '{print $8}'`
    echo "* hotspot.py (GCAP) version: $var"; versions=`echo "$versions, \"hotspot.py (GCAP)\": \"$var\""`
    var=`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* samtools version: $var"; versions=`echo "$versions, \"samtools\": \"$var\""`
    var=`bedops --version 2>&1 | grep version | awk '{print $2}'`
    echo "* bedops version: $var"; versions=`echo "$versions, \"bedops\": \"$var\""`
    var=`bedmap --version 2>&1 | grep version | awk '{print $2}'`
    echo "* bedmap (bedops) version: $var"; versions=`echo "$versions, \"bedmap (bedops)\": \"$var\""`
    var=`sort-bed --version 2>&1 | grep version | awk '{print $2}'`
    echo "* sort-bed (bedops version: $var"; versions=`echo "$versions, \"sort-bed (bedops\": \"$var\""`
    var=`starch --version 2>&1 | grep version | awk '{print $3}'`
    echo "* starch (bedops) version: $var"; versions=`echo "$versions, \"starch (bedops)\": \"$var\""`
    var=`starchcat --version 2>&1 | grep version | awk '{print $3}'`
    echo "* starchcat (bedops) version: $var"; versions=`echo "$versions, \"starchcat (bedops)\": \"$var\""`
    var=`unstarch --version 2>&1 | grep version | awk '{print $3}'`
    echo "* unstarch (bedops) version: $var"; versions=`echo "$versions, \"unstarch (bedops)\": \"$var\""`
    var=`bedtools --version 2>&1 | awk '{print $2}'`
    echo "* bedtools version: $var"; versions=`echo "$versions, \"bedtools\": \"$var\""`
    var=`bamToBed -h 2>&1 | grep Version | awk '{print $2}'`
    echo "* bamToBed (bedtools) version: $var"; versions=`echo "$versions, \"bamToBed (bedtools)\": \"$var\""`
    var=`shuffleBed -h 2>&1 | grep Version | awk '{print $2}'`
    echo "* shuffleBed (bedtools) version: $var"; versions=`echo "$versions, \"shuffleBed (bedtools)\": \"$var\""`
    versions=`echo $versions }`
    echo "*****"

    echo "* Value of bam_no_chrM:   '$bam_no_chrM'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of read_length: '$read_length'"
    echo "* Value of genome: '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_no_chrM" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_no_chrM" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    sample_root="${bam_root}_5M_sample"
    echo "* out: sample file: '${sample_root}.bam'"

    dx download "$chrom_sizes" -o chromSizes.txt
    # sort-bed is important!
    grep -v chrM chromSizes.txt | awk '{printf "%s\t0\t%s\t%s\n",$1,$2,$1}' | sort-bed - > ${genome}.chromInfo.bed

    read_size=`dx describe "$bam_no_chrM" --details --json | grep "\"average length:\"" | awk '{print $3}' | tr -d ,`
    echo "* Found read size: '$read_size'"
    if [ "$read_size" != "" ] && [ "$read_length" -ne "$read_size" ]; then
        echo "* WARNING Read length ($read_length) does not match discovered read size ($read_size)."
    fi

    # TODO: Need to make bam.bai?
    echo "* Indexing bam..."
    set -x
    samtools index ${bam_root}.bam
    set +x
    
    echo "* Sampling bam..."
    set -x
    edwBamStats -sampleBamSize=5000000 -u4mSize=5000000 -sampleBam=${sample_root}.bam ${bam_root}.bam ${sample_root}_stats.txt
    cat ${sample_root}_stats.txt
    samtools index ${sample_root}.bam
    set +x
    
    echo "* Running hotspot.py..."
    set -x
    mkdir tmp
    mkdir out
    cp ${genome}.chromInfo.bed hotspot/hotspot-distr/data/
    # May also need to do something about "Satellite.${genome}.bed"
    #cp /usr/bin/Satellite.${genome}.bed hotspot/hotspot-distr/data   # hg19 version already there!
    mappable=${genome}.K${read_length}.mappable_only
    wget http://www.uwencode.org/proj/hotspot/${mappable}.bed -O hotspot/hotspot-distr/data/${mappable}.bed >> install.log 2>&1
    python2.7 /usr/bin/hotspot.py -o hotspot/hotspot-distr/ ${sample_root}.bam $genome DNase-seq $read_length tmp out
    # hanging on: /home/dnanexus/hotspot/hotspot-distr/pipeline-scripts/run_generate_random_lib

    cp tmp/${sample_root}.spot.out ${sample_root}_hotspot_qc.txt
    set +x

    echo "* Prepare metadata..."
    meta=`echo \"edwBamStats\": { `
    # alignedBy bwa
    var=`grep "^alignedBy" ${sample_root}_stats.txt | awk '{printf "\"alignedBy\": \"%s\"", $2}'`
    meta=`echo $meta $var`
    # isPaired 0
    var=`grep "^isPaired" ${sample_root}_stats.txt | awk '{printf "\"isPaired\": %d", $2}'`
    meta=`echo $meta, $var`
    # isSortedByTarget 1
    var=`grep "^isSortedByTarget" ${sample_root}_stats.txt | awk '{printf "\"isSortedByTarget\": %d", $2}'`
    meta=`echo $meta, $var`
    # readCount 2286836
    var=`grep "^readCount" ${sample_root}_stats.txt | awk '{printf "\"readCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # readBaseCount 82326096
    var=`grep "^readBaseCount" ${sample_root}_stats.txt | awk '{printf "\"readBaseCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # mappedCount 2286836
    var=`grep "^mappedCount" ${sample_root}_stats.txt | awk '{printf "\"mappedCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # uniqueMappedCount 2285489
    var=`grep "^uniqueMappedCount" ${sample_root}_stats.txt | awk '{printf "\"uniqueMappedCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeMean 36
    var=`grep "^readSizeMean" ${sample_root}_stats.txt | awk '{printf "\"readSizeMean\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeStd 0
    var=`grep "^readSizeStd" ${sample_root}_stats.txt | awk '{printf "\"readSizeStd\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeMin 36
    var=`grep "^readSizeMin" ${sample_root}_stats.txt | awk '{printf "\"readSizeMin\": %d", $2}'`
    meta=`echo $meta, $var`
    # readSizeMax 36
    var=`grep "^readSizeMax" ${sample_root}_stats.txt | awk '{printf "\"readSizeMax\": %d", $2}'`
    meta=`echo $meta, $var`
    # u4mReadCount 2285489
    var=`grep "^u4mReadCount" ${sample_root}_stats.txt | awk '{printf "\"u4mReadCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # u4mUniquePos 2213397
    var=`grep "^u4mUniquePos" ${sample_root}_stats.txt | awk '{printf "\"u4mUniquePos\": %d", $2}'`
    meta=`echo $meta, $var`
    # u4mUniqueRatio 0.968457
    var=`grep "^u4mUniqueRatio" ${sample_root}_stats.txt | awk '{printf "\"u4mUniqueRatio\": %s", $2}'`
    meta=`echo $meta, $var`
    # targetSeqCount 25
    var=`grep "^targetSeqCount" ${sample_root}_stats.txt | awk '{printf "\"targetSeqCount\": %d", $2}'`
    meta=`echo $meta, $var`
    # targetBaseCount 3095693983
    var=`grep "^targetBaseCount" ${sample_root}_stats.txt | awk '{printf "\"targetBaseCount\": %d", $2}'`
    meta=`echo $meta, $var`
    meta=`echo $meta }`
    qc_sampled=$meta

    meta=`echo \"hotspot\": { `
    #total tags  hotspot tags    SPOT
    # 2255195       1083552  0.4804
    var=`tail -1 ${sample_root}_hotspot_qc.txt | awk '{printf "\"total tags\": %d, \"hotspot tags\": %d, \"SPOT\": %d", $1,$2,$3}'`
    meta=`echo $meta $var`
    meta=`echo $meta }`
    qc_sampled=`echo $qc_sampled, $meta`

    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bam_5M_sample=$(dx upload ${sample_root}.bam --details "{ $qc_sampled }" --property QC="$qc_sampled" --property SW="$versions" --brief)
    #bam_5M_sample=$(dx upload ${sample_root}.bam --brief)
    hotspot_5M_sample_qc=$(dx upload ${sample_root}_hotspot_qc.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bam_5M_sample "$bam_5M_sample" --class=file
    dx-jobutil-add-output hotspot_5M_sample_qc "$hotspot_5M_sample_qc" --class=file
    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
