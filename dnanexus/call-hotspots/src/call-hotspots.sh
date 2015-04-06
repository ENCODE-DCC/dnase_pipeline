#!/bin/bash
# call-hotspots 0.0.1

main() {
    echo "* Installing hotspot and dependencies (gsl)..." 2>&1 | tee -a install.log
    exe_dir="`pwd`"
    echo "> ls -l /usr/local/lib" 2>&1 | tee -a install.log
    set -x
    ls -l /usr/local/lib 2>&1 | tee -a install.log
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz -O gsl.tgz 2>&1 | tee -a install.log
    mkdir gsl
    tar -xzf gsl.tgz -C gsl --strip-components=1
    cd gsl
    ./configure 2>&1 | tee -a install.log
    make 2>&1 | tee -a install.log
    sudo make install 2>&1 | tee -a install.log
    set +x; echo "> ls -l /usr/local/lib" 2>&1 | tee -a install.log; set -x
    ls -l /usr/local/lib 2>&1 | tee -a install.log
    set +x; echo "gsl-config --libs" 2>&1 | tee -a install.log; set -x
    gsl-config --libs 2>&1 | tee -a install.log
    set +x; echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}" 2>&1 | tee -a install.log; set -x
    export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
    set +x; echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}" 2>&1 | tee -a install.log; set -x
    cd ..
    wget https://github.com/rthurman/hotspot/archive/v4.1.0.tar.gz -O hotspot.tgz 2>&1 | tee -a install.log
    mkdir hotspot
    tar -xzf hotspot.tgz -C hotspot --strip-components=1
    cd hotspot/hotspot-distr/hotspot-deploy
    make 2>&1 | tee -a install.log
    set +x; echo "> ls -l bin" 2>&1 | tee -a install.log; set -x
    ls -l bin 2>&1 | tee -a install.log
    # Can either put bin in path or copy contents of bin to /usr/bin
    export PATH=${exe_dir}/hotspot/hotspot-distr/hotspot-deploy/bin:${PATH}
    #cp bin/* /usr/bin
    set +x; echo "PATH: $PATH" 2>&1 | tee -a install.log
            echo "* === hotspot ===" 2>&1 | tee -a install.log; set -x
    bin/hotspot 2>&1 | tee -a install.log
    set +x; echo "* ===============" 2>&1 | tee -a install.log; set -x
    hotspot 2>&1 | tee -a install.log
    set +x; echo "* ===============" 2>&1 | tee -a install.log; set -x
    cd ../../../
    # additional executables in resources/usr/bin
    set +x
    
    echo "*****"
    echo "* Running: call-hotspots.sh v0.0.1"
    #overkill: echo "* gsl: "`gsl/configure --version 2> /dev/null | grep gsl | awk '{print $3}'` | tee -a install.log
    echo "* hotspot: "`hotspot 2>&1 | grep HotSpot | awk '{print $1}'`
    echo "* bedops: "`bedops --version 2>&1 | grep version | awk '{print $2}'`
    #echo "* bedmap (bedops): "`bedmap --version 2>&1 | grep version | awk '{print $2}'`
    #echo "* sort-bed (bedops): "`sort-bed --version 2>&1 | grep version | awk '{print $2}'`
    #echo "* starch (bedops): "`starch --version 2>&1 | grep version | awk '{print $3}'`
    #echo "* starchcat (bedops): "`starchcat --version 2>&1 | grep version | awk '{print $3}'`
    echo "* unstarch (bedops): "`unstarch --version 2>&1 | grep version | awk '{print $3}'`
    #not installed echo "* GCAP: (unversioned) "`GCAP/gcap/GCAP -h 2> /dev/null | grep Global`
    echo "* hotspot.py (GCAP): "`python2.7 hotspot.py -h | grep Version | awk '{print $8}'`
    echo "* bedToBigBed: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{print $2$3}'`
    echo "* bedGraphToWig: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    echo "* bedGraphPack: "`./bedGraphPack 2>&1 | grep "bedGraphPack v" | awk '{print $2}'`
    echo "* intersectBed (bedtools): "`intersectBed 2>&1 | grep Version | awk '{print $2}'`
    echo "* samtools: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "*****"

    echo "* Value of bam_input    '$bam_input'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of read_length: '$read_length'"
    echo "* Value of genome:      '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_input" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_input" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    narrowPeak_root="${bam_root}_narrowPeak_hotspot"
    broadPeak_root="${bam_root}_broadPeak_hotspot"
    signal_root="${bam_root}_signal_hotspot"
    echo "* out: narrowPeak files: '${narrowPeak_root}.bed' and '${narrowPeak_root}.bb'"
    echo "* out: broadPeak files: '${broadPeak_root}.bed' and '${broadPeak_root}.bb'"
    echo "* out: signal file: '${signal_root_root}.bw'"

    dx download "$chrom_sizes" -o chromSizes.txt
    grep -v chrM chromSizes.txt | sort | awk '{printf "%s\t0\t%s\t%s\n",$1,$2,$1}' > ${genome}.chromInfo.bed

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
    cp /usr/bin/${genome}.K${read_length}.mappable_only.bed hotspot/hotspot-distr/data
    cp /usr/bin/${genome}.K${read_length}.mappable_only.starch hotspot/hotspot-distr/data
    #cp /usr/bin/Satellite.${genome}.bed hotspot/hotspot-distr/data   # hg19 version already there!
    python2.7 hotspot.py hotspot/hotspot-distr/ ${bam_root}.bam $genome DNase-seq $read_length tmp out
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



    # TODO: Collect QC from: ?
    
    echo "* Prepare metadata..."
    meta=`echo \"hotspot_out\": { `
    ## # 2142994 + 0 in total (QC-passed reads + QC-failed reads)
    ## var=`grep "QC-passed reads" ${bam_filtered_root}_qc.txt | awk '{printf "\"total\": %d, \"total_qc_failed\": %d", $1,$3}'`
    ## meta=`echo $meta $var`
    ## # 0 + 0 duplicates
    ## var=`grep -w duplicates ${bam_filtered_root}_qc.txt | awk '{printf "\"duplicates\": %d, \"duplicates_qc_failed\": %d", $1,$3}'`
    ## meta=`echo $meta, $var`
    ## # 2046212 + 0 mapped (95.48%:-nan%)
    ## var=`grep -w mapped ${bam_filtered_root}_qc.txt | awk '{printf "\"mapped\": %d, \"mapped_qc_failed\": %d", $1,$3}'`
    ## meta=`echo $meta, $var`
    ## var=`grep -w mapped ${bam_filtered_root}_qc.txt | awk '{print $5}' | tr ":" " " | awk '{print $1}' | tr -d "("`
    ## meta=`echo $meta, \"mapped_pct\": \"$var\"`
    meta=`echo $meta }`



    echo "* Upload results..."
    # NOTE: adding meta 'details' ensures json is valid.  But details are not updatable so rely on QC property
    #details=`echo { $meta }`
    #bam_filtered=$(dx upload ${bam_filtered_root}.bam --details "$details" --property QC="$meta" --brief)
    bed_hotspot_narrowPeak=$(dx upload ${narrowPeak_root}.bed --brief)
    bed_hotspot_broadPeak=$(dx upload ${broadPeak_root}.bed --brief)
    bb_hotspot_narrowPeak=$(dx upload ${narrowPeak_root}.bb --brief)
    bb_hotspot_broadPeak=$(dx upload ${broadPeak_root}.bb --brief)
    bw_hotspot_signal=$(dx upload ${signal_root}.bw --brief)
    bam_hotspot_qc=$(dx upload ${bam_root}_hotspot_qc.txt --brief)

    dx-jobutil-add-output bed_hotspot_narrowPeak "$bed_hotspot_narrowPeak" --class=file
    dx-jobutil-add-output bed_hotspot_broadPeak "$bed_hotspot_broadPeak" --class=file
    dx-jobutil-add-output bb_hotspot_narrowPeak "$bb_hotspot_narrowPeak" --class=file
    dx-jobutil-add-output bb_hotspot_broadPeak "$bb_hotspot_broadPeak" --class=file
    dx-jobutil-add-output bw_hotspot_signal "$bw_hotspot_signal" --class=file
    dx-jobutil-add-output bam_hotspot_qc "$bam_hotspot_qc" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
