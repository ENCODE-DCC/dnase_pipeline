#!/bin/bash
# sample-hotspots 0.0.1

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
    echo "* Running: sample-hotspots.sh v0.0.1"
    #overkill: echo "* gsl version: "`gsl/configure --version 2> /dev/null | grep gsl | awk '{print $3}'` | tee -a install.log
    echo "* hotspot version: "`hotspot 2>&1 | grep HotSpot | awk '{print $1}'`
    echo "* bedops version: "`bedops --version 2>&1 | grep version | awk '{print $2}'`
    #echo "* bedmap (bedops) version: "`bedmap --version 2>&1 | grep version | awk '{print $2}'`
    #echo "* sort-bed (bedops) version: "`sort-bed --version 2>&1 | grep version | awk '{print $2}'`
    #echo "* starch (bedops) version: "`starch --version 2>&1 | grep version | awk '{print $3}'`
    #echo "* starchcat (bedops) version: "`starchcat --version 2>&1 | grep version | awk '{print $3}'`
    #echo "* unstarch (bedops) version: "`unstarch --version 2>&1 | grep version | awk '{print $3}'`
    #not installed echo "* GCAP version: (unversioned) "`GCAP/gcap/GCAP -h 2> /dev/null | grep Global`
    echo "* hotspot.py (GCAP) version: "`python2.7 hotspot.py -h | grep Version | awk '{print $8}'`
    #echo "* bedToBigBed version: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{print $2$3}'`
    #echo "* bedGraphToWig version: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    #echo "* bedGraphPack version: "`./bedGraphPack 2>&1 | grep "bedGraphPack v" | awk '{print $2}'`
    echo "* intersectBed (bedtools) version: "`intersectBed 2>&1 | grep Version | awk '{print $2}'`

    echo "* samtools version: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* edwBamFilter version: "`edwBamFilter 2>&1 | grep "edwBamFilter v" | awk '{print $2}'`
    echo "* edwBamStats version: "`edwBamStats 2>&1 | grep "edwBamStats v" | awk '{print $2}'`
    echo "* phantompeakqualtools version: "`grep Version phantompeakqualtools/README.txt | awk '{print $2}'`
    echo "* bedtools version: "`bedtools --version 2>&1 | awk '{print $2}'`
    echo "*****"

    echo "* Value of bam_input:   '$bam_input'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of read_length: '$read_length'"
    echo "* Value of genome: '$genome'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_input" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_input" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    sample_root="${bam_root}_5M_sample"
    echo "* out: sample file: '${sample_root}.bam'"

    dx download "$chrom_sizes" -o chromSizes.txt
    grep -v chrM chromSizes.txt | sort | awk '{printf "%s\t0\t%s\t%s\n",$1,$2,$1}' > ${genome}.chromInfo.bed

    # TODO: Need to make bam.bai?
    echo "* Indexing bam..."
    set -x
    samtools index ${bam_root}.bam
    set +x
    
    echo "* Sampling bam..."
    set -x
    edwBamStats -sampleBamSize=5000000 -u4mSize=5000000 -sampleBam=${sample_root}.bam ${bam_root}.bam tmp.ra
    samtools index ${sample_root}.bam
    set +x

    echo "* Running hotspot.py..."
    set -x
    mkdir tmp
    mkdir out
    cp ${genome}.chromInfo.bed hotspot/hotspot-distr/data/
    cp /usr/bin/${genome}.K${read_length}.mappable_only.bed hotspot/hotspot-distr/data
    cp /usr/bin/${genome}.K${read_length}.mappable_only.starch hotspot/hotspot-distr/data
    #cp /usr/bin/Satellite.${genome}.bed hotspot/hotspot-distr/data   # hg19 version already there!
    python2.7 hotspot.py -o hotspot/hotspot-distr/ ${sample_root}.bam $genome DNase-seq $read_length tmp out
    cp tmp/${sample_root}.spot.out ${sample_root}_hotspot_qc.txt
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
    bam_5M_sample=$(dx upload ${sample_root}.bam --brief)
    hotspot_5M_sample_qc=$(dx upload ${sample_root}_hotspot_qc.txt --brief)

    dx-jobutil-add-output bam_5M_sample "$bam_5M_sample" --class=file
    dx-jobutil-add-output hotspot_5M_sample_qc "$hotspot_5M_sample_qc" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
