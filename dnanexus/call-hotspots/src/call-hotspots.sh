#!/bin/bash
# call-hotspots 0.1.0

script_name="call-hotspots.sh"
script_ver="0.1.0"
#script_name=`head -2 $0 | tail -1 | awk '{ print $2 }'`
#script_ver=`head -2 $0 | tail -1 | awk '{ print $3 }'`

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
    var=`intersectBed 2>&1 | grep Version | awk '{print $2}'`
    echo "* intersectBed (bedtools) version: $var"; versions=`echo "$versions, \"intersectBed (bedtools)\": \"$var\""`
    var=`shuffleBed -h 2>&1 | grep Version | awk '{print $2}'`
    echo "* shuffleBed (bedtools) version: $var"; versions=`echo "$versions, \"shuffleBed (bedtools)\": \"$var\""`
    var=`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{print $2$3}'`
    echo "* bedToBigBed version: $var"; versions=`echo "$versions, \"bedToBigBed\": \"$var\""`
    var=`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    echo "* bedGraphToBigWig version: $var"; versions=`echo "$versions, \"bedGraphToBigWig\": \"$var\""`    
    var=`bedGraphPack 2>&1 | grep "bedGraphPack v" | awk '{print $2}'`
    echo "* bedGraphPack version: $var"; versions=`echo "$versions, \"bedGraphPack\": \"$var\""`
    versions=`echo $versions }`
    echo "*****"

    echo "* Value of bam_no_chrM: '$bam_no_chrM'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of read_length: '$read_length'"
    echo "* Value of genome:      '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_no_chrM" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_no_chrM" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    narrowPeak_root="${bam_root}_narrowPeak_hotspot"
    broadPeak_root="${bam_root}_broadPeak_hotspot"
    signal_root="${bam_root}_signal_hotspot"
    echo "* out: narrowPeak files: '${narrowPeak_root}.bed' and '${narrowPeak_root}.bb'"
    echo "* out: broadPeak files: '${broadPeak_root}.bed' and '${broadPeak_root}.bb'"
    echo "* out: signal file: '${signal_root_root}.bw'"

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
        
    echo "* Running hotspot.py..."
    set -x
    mkdir tmp
    mkdir out
    cp ${genome}.chromInfo.bed hotspot/hotspot-distr/data/
    # May also need to do something about "Satellite.${genome}.bed"
    #cp /usr/bin/Satellite.${genome}.bed hotspot/hotspot-distr/data   # hg19 version already there!
    mappable=${genome}.K${read_length}.mappable_only
    wget http://www.uwencode.org/proj/hotspot/${mappable}.bed -O hotspot/hotspot-distr/data/${mappable}.bed >> install.log 2>&1
    python2.7 /usr/bin/hotspot.py hotspot/hotspot-distr/ ${bam_root}.bam $genome DNase-seq $read_length tmp out
    cp tmp/${bam_root}.spot.out ${bam_root}_hotspot_qc.txt
    set +x

    echo "* Generating narrowPeaks..."
    set -x
    # Convert output into ENCODE formats,
    # Convert narrowPeaks.bed and several stray columns to narrowPeaks.bigBed in $4
    paste out/narrowPeaks.bed out/narrowPeaks.dens out/narrowPeaks.pval | \
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "p" NR, 0, ".", $4, $5, -1, -1}' - | \
        sort -k1,1 -k2,2n - > ${narrowPeak_root}.bed
    bedToBigBed -as=/usr/bin/narrowPeak.as -type=bed6+4 ${narrowPeak_root}.bed chromSizes.txt ${narrowPeak_root}.bb
    wc -l ${narrowPeak_root}.bed > ${narrowPeak_root}_qc.txt
    set +x

    echo "* Generating broadPeaks..."
    set -x
    # Convert broadPeaks.bed and broadPeaks.pVal to broadPeaks.bigBed in $5
    paste out/broadPeaks.bed out/broadPeaks.pval | \
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "hot" NR, 0, ".", $5, $6, -1}' - | \
        sort -k1,1 -k2,2n - > ${broadPeak_root}.bed
    bedToBigBed -as=/usr/bin/broadPeak.as -type=bed6+3 ${broadPeak_root}.bed chromSizes.txt ${broadPeak_root}.bb
    wc -l ${broadPeak_root}.bed > ${broadPeak_root}_qc.txt
    set +x

    echo "* Generating signal..."
    set -x
    # Convert starched bedGraph to mappable-only bigWig in $6
    unstarch out/density.bed.starch > tmp/tmp.bed
    intersectBed -a tmp/tmp.bed -b hotspot/hotspot-distr/data/{$genome}.K${read_length}.mappable_only.bed -f 1.00 | \
        cut -f 1,2,3,5 | bedGraphPack stdin ${signal_root_root}.bedGraph
    bedGraphToBigWig ${signal_root_root}.bedGraph chromSizes.txt ${signal_root_root}.bw
    set +x

    echo "* Prepare metadata..."
    meta=`echo \"hotspot\": { `
    #total tags  hotspot tags    SPOT
    # 2255195       1083552  0.4804
    var=`tail -1 ${bam_root}_hotspot_qc.txt | awk '{printf "\"total tags\": %d, \"hotspot tags\": %d, \"SPOT\": %d", $1,$2,$3}'`
    meta=`echo $meta $var`
    meta=`echo $meta }`
    qc_hotspot=$meta
    
    meta=`echo \"peak_counts\": { `
    var=`cat ${narrowPeak_root}_qc.txt | awk '{printf "\"hotspot count\": %d", $1}'`
    qc_spots=`echo $qc_hotspot, $var`
    meta=`echo $meta $var`
    var=`cat ${broadPeak_root}_qc.txt | awk '{printf "\"regions count\": %d", $1}'`
    qc_regions=`echo $qc_hotspot, $var`
    meta=`echo $meta, $var`
    meta=`echo $meta }`
    qc_hotspot=`echo $qc_hotspot, $meta`
    
    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    bed_hotspot_narrowPeak=$(dx upload ${narrowPeak_root}.bed --details "{ $qc_spots }"   --property QC="$qc_spots"   --property SW="$versions" --brief)
    bed_hotspot_broadPeak=$(dx upload ${broadPeak_root}.bed   --details "{ $qc_regions }" --property QC="$qc_regions" --property SW="$versions" --brief)
    bb_hotspot_narrowPeak=$(dx upload ${narrowPeak_root}.bb   --details "{ $qc_spots }"   --property QC="$qc_spots"   --property SW="$versions" --brief)
    bb_hotspot_broadPeak=$(dx upload ${broadPeak_root}.bb     --details "{ $qc_regions }" --property QC="$qc_regions" --property SW="$versions" --brief)
    bw_hotspot_signal=$(dx upload ${signal_root}.bw           --details "{ $qc_hotspot }" --property QC="$qc_hotspot" --property SW="$versions" --brief)
    bam_hotspot_qc=$(dx upload ${bam_root}_hotspot_qc.txt --property SW="$versions" --brief)

    dx-jobutil-add-output bed_hotspot_narrowPeak "$bed_hotspot_narrowPeak" --class=file
    dx-jobutil-add-output bed_hotspot_broadPeak "$bed_hotspot_broadPeak" --class=file
    dx-jobutil-add-output bb_hotspot_narrowPeak "$bb_hotspot_narrowPeak" --class=file
    dx-jobutil-add-output bb_hotspot_broadPeak "$bb_hotspot_broadPeak" --class=file
    dx-jobutil-add-output bw_hotspot_signal "$bw_hotspot_signal" --class=file
    dx-jobutil-add-output bam_hotspot_qc "$bam_hotspot_qc" --class=file
    dx-jobutil-add-output metadata "$versions" --class=string

    echo "* Finished."
}
