#!/bin/bash
# dnase-install 0.0.1

main() {
    # Some executables in resources/usr/bin
    echo "*****"
    echo "* Running: dnase-install.sh v0.0.1" | tee -a install.log
    echo "*****"
    exe_dir="`pwd`"
    gentoken="${genome}_${gender}"
    
    echo "* Installing python3 python3-dev..." 2>&1 | tee -a install.log
    wget https://www.python.org/ftp/python/3.3.3/Python-3.3.3.tgz -O python3.tgz 2>&1 | tee -a install.log
    mkdir python3
    tar -xzf python3.tgz -C python3 --strip-components=1
    cd python3
    ./configure 2>&1 | tee -a install.log
    make 2>&1 | tee -a install.log
    #make test 2>&1 | tee -a install.log
    sudo make install 2>&1 | tee -a install.log
    alias python3=python3.3
    cd ..
    #sudo apt-get --yes install python3 python3-dev 2>&1 | tee -a install.log

    echo "* Installing setuptools..." 2>&1 | tee -a install.log
    wget https://bootstrap.pypa.io/ez_setup.py -O - | sudo python3 2>&1 | tee -a install.log

    #echo "* Installing R 3.1..." 2>&1 | tee -a install.log
    ##sudo apt-get update 2>&1 | tee -a install.log
    ##sudo apt-get install r-base 2>&1 | tee -a install.log
    ##sudo apt-get update 2>&1 | tee -a install.log
    ##sudo apt-get upgrade 2>&1 | tee -a install.log   ### Danger screw up with interactive msql upgrade
    #wget http://cran.cnr.berkeley.edu/src/base/R-3/R-3.1.1.tar.gz -O r3.tgz 2>&1 | tee -a install.log
    #mkdir r3
    #tar -xzf r3.tgz -C r3 --strip-components=1
    #cd r3
    #./configure --with-x=no 2>&1 | tee -a install.log
    #JAVA_HOME=`whereis java | awk '{print $2}'`
    #make 2>&1 | tee -a install.log
    ##make check 2>&1 | tee -a install.log
    #sudo make install 2>&1 | tee -a install.log
    #cd ..

    echo "* Installing phantompeakqualtools and dependencies (caTools, snow)..." 2>&1 | tee -a install.log
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

    echo "* Installing hotspot and dependencies (gsl)..." 2>&1 | tee -a install.log
    echo "> ls -l /usr/local/lib" 2>&1 | tee -a install.log
    ls -l /usr/local/lib 2>&1 | tee -a install.log
    wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz -O gsl.tgz 2>&1 | tee -a install.log
    mkdir gsl
    tar -xzf gsl.tgz -C gsl --strip-components=1
    cd gsl
    ./configure 2>&1 | tee -a install.log
    make 2>&1 | tee -a install.log
    sudo make install 2>&1 | tee -a install.log
    echo "> ls -l /usr/local/lib" 2>&1 | tee -a install.log
    ls -l /usr/local/lib 2>&1 | tee -a install.log
    echo "gsl-config --libs" 2>&1 | tee -a install.log
    gsl-config --libs 2>&1 | tee -a install.log
    echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}" 2>&1 | tee -a install.log
    export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}
    echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}" 2>&1 | tee -a install.log
    cd ..
    wget https://github.com/rthurman/hotspot/archive/v4.1.0.tar.gz -O hotspot.tgz 2>&1 | tee -a install.log
    mkdir hotspot
    tar -xzf hotspot.tgz -C hotspot --strip-components=1
    cd hotspot/hotspot-distr/hotspot-deploy
    make 2>&1 | tee -a install.log
    ls -l bin 2>&1 | tee -a install.log
    # Can either put bin in path or copy contents of bin to /usr/bin
    export PATH=${PATH}:${exe_dir}/hotspot/hotspot-distr/hotspot-deploy/bin
    #cp bin/* /usr/bin
    echo "PATH: $PATH" 2>&1 | tee -a install.log
    echo "* === hotspot ===" 2>&1 | tee -a install.log
    bin/hotspot 2>&1 | tee -a install.log
    echo "* ===============" 2>&1 | tee -a install.log
    hotspot 2>&1 | tee -a install.log
    echo "* ===============" 2>&1 | tee -a install.log
    echo "pwd "`pwd` 2>&1 | tee -a install.log
    cd ../../../
    echo "pwd "`pwd` 2>&1 | tee -a install.log

    echo "* Installing macs2 and dependencies (numpy)..." 2>&1 | tee -a install.log
    wget http://sourceforge.net/projects/numpy/files/NumPy/1.6.2/numpy-1.6.2.tar.gz/download -O numpy.tgz 2>&1 | tee -a install.log
    mkdir numpy
    tar -xzf numpy.tgz -C numpy --strip-components=1
    cd numpy
    python2.7 setup.py build --fcompiler=gnu95 2>&1 | tee -a install.log
    sudo python2.7 setup.py install 2>&1 | tee -a install.log
    cd ..
    wget https://pypi.python.org/packages/source/M/MACS2/MACS2-2.0.10.20131216.tar.gz -O macs2.tgz 2>&1 | tee -a install.log
    mkdir macs2
    tar -xzf macs2.tgz -C macs2 --strip-components=1
    cd macs2
    sudo python2.7 setup.py install 2>&1 | tee -a install.log
    cd ..

    echo "* Installing GCAP..." 2>&1 | tee -a install.log
    git clone https://github.com/qinqian/GCAP.git
    #wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz -O gcap.tgz 2>&1 | tee -a install.log
    #mkdir GCAP
    #tar -xzf gcap.tgz -C GCAP --strip-components=1
    cd GCAP
    # edit gcap.conf
    ##mv -f gcap.conf gcap.conf.orig
    echo "[tool]" > gcap.conf
    echo "peak_calling = ${exe_dir}/hotspot/hotspot-distr" >> gcap.conf
    echo "## for strand correlation NSC, RSA" >> gcap.conf
    echo "spp = ${exe_dir}/phantompeakqualtools/run_spp.R" >> gcap.conf
    echo " " >> gcap.conf
    echo "[${gentoken}]" >> gcap.conf
    echo "genome_index = ${exe_dir}/${gentoken}.fa" >> gcap.conf
    echo "chrom_len = ${exe_dir}/${gentoken}_chrom.sizes" >> gcap.conf
    echo "chrom_bed = ${exe_dir}/${gentoken}_chrom.bed" >> gcap.conf
    echo "* === gcap.conf ===" 2>&1 | tee -a install.log
    cat gcap.conf 2>&1 | tee -a install.log
    echo "* =================" 2>&1 | tee -a install.log
    sudo python3 setup.py install 2>&1 | tee -a install.log
    cd ..
    #echo "* === GCAP -h ===" 2>&1 | tee -a install.log
    #GCAP/gcap/GCAP -h 2>&1 | tee -a install.log
    #echo "* =====================" 2>&1 | tee -a install.log
  
    echo "*****"
    echo "* Running: dnase-install.sh v0.0.1"
    echo "* python3: "`python3 --version 2>&1 | awk '{print $2}'`
    echo "* python2.7: "`python2.7 --version 2>&1 | awk '{print $2}'`
    echo "* R: "`R --version | grep "R version" | awk '{print $3,$4}'`
    echo "* phantompeakqualtools: "`grep Version phantompeakqualtools/README.txt | awk '{print $2}'`
    VER=`grep caTools_ phantompeakqualtools/installPkgs.R | tr \_ " " | awk '{print $2}'`
    echo "* caTools: ${VER%.tar*}"
    VER=`grep snow_ phantompeakqualtools/installPkgs.R | tr \_ " " | awk '{print $2}'`
    echo "* snow: ${VER%.tar*}"
    echo "* bwa: "`bwa 2>&1 | grep Version | awk '{print $2}'`
    echo "* fastqStatsAndSubsample: "`fastqStatsAndSubsample 2>&1 | grep "fastqStatsAndSubsample v" | awk '{print $2}'`
    echo "* edwBamFilter: "`edwBamFilter 2>&1 | grep "edwBamFilter v" | awk '{print $2}'`
    echo "* edwBamStats: "`edwBamStats 2>&1 | grep "edwBamStats v" | awk '{print $2}'`
    echo "* edwComparePeaks: (unversioned) "`edwComparePeaks 2>&1 | grep "edwComparePeaks -" | awk '{print $3,$4,$5,$6}'`
    echo "* bigWigCorrelate: (unversioned) "`bigWigCorrelate 2>&1 | grep "bigWigCorrelate -" | awk '{print $3,$4,$5}'`
    echo "* bedToBigBed: "`bedToBigBed 2>&1 | grep "bedToBigBed v" | awk '{print $2$3}'`
    echo "* bigBedToBed: "`bigBedToBed 2>&1 | grep "bigBedToBed v" | awk '{print $2}'`
    echo "* bedGraphToWig: "`bedGraphToBigWig 2>&1 | grep "bedGraphToBigWig v" | awk '{print $2$3}'`
    echo "* bedGraphPack: "`./bedGraphPack 2>&1 | grep "bedGraphPack v" | awk '{print $2}'`
    echo "* samtools: "`samtools 2>&1 | grep Version | awk '{print $2}'`
    echo "* bedops: "`bedops --version 2>&1 | grep version | awk '{print $2}'`
    echo "* bedmap (bedops): "`bedmap --version 2>&1 | grep version | awk '{print $2}'`
    echo "* sort-bed (bedops): "`sort-bed --version 2>&1 | grep version | awk '{print $2}'`
    echo "* starch (bedops): "`starch --version 2>&1 | grep version | awk '{print $3}'`
    echo "* starchcat (bedops): "`starchcat --version 2>&1 | grep version | awk '{print $3}'`
    echo "* unstarch (bedops): "`unstarch --version 2>&1 | grep version | awk '{print $3}'`
    echo "* gsl: "`gsl/configure --version 2> /dev/null | grep gsl | awk '{print $3}'`
    echo "* hotspot: "`hotspot 2>&1 | grep HotSpot | awk '{print $1}'`
    echo "* numpy: "`grep Version numpy/PKG-INFO | grep -v Metadata | awk '{print $2}'`
    echo "* macs2: "`macs2 --version 2>&1 | awk '{print $2" "$3}'`
    echo "* bedtools: "`bedtools --version 2>&1 | awk '{print $2}'`
    echo "* bedToBam (bedtools): "`bedToBam 2>&1 | grep Version | awk '{print $2}'`
    echo "* intersectBed (bedtools): "`intersectBed 2>&1 | grep Version | awk '{print $2}'`
    echo "* GCAP: (unversioned) "`GCAP/gcap/GCAP -h 2> /dev/null | grep Global`
    echo "* hotspot.py (GCAP): "`python hotspot.py -h | grep Version | awk '{print $8}'`
    echo "*****"
    # from GCAP readme.md:
    #python(for pipeline)	3.3	https://www.python.org/ftp/python/3.3.3/Python-3.3.3.tgz  INSTALL
    #python(for hotspot)	2.7	http://www.python.org/ftp/python/2.7.6/Python-2.7.6.tgz   EXPECT
    #R	3.1.1	http://cran.cnr.berkeley.edu/src/base/R-3/R-3.1.1.tar.gz                  EXPECT
    #bwa	0.7.7	https://github.com/lh3/bwa                                            COPIED TO resources/usr/bin (bwa)
    #fastqStatsAndSubsample	2	https://github.com/ENCODE-DCC/kentUtils                   COPIED TO resources/usr/bin
    #edwBamFilter	2	https://github.com/ENCODE-DCC/kentUtils                           COPIED TO resources/usr/bin
    #edwBamStats	2	https://github.com/ENCODE-DCC/kentUtils                           COPIED TO resources/usr/bin
    #samtools	0.2.0	https://github.com/samtools/samtools                              COPIED TO resources/usr/bin
    #bedops	2.4.2	https://github.com/bedops/bedops/releases/download/v2.4.2/bedops_linux_x86_64-v2.4.2.tar.bz2   COPIED TO resources/usr/bin  (bedops, bedmap, sort-bed, starch, starchcat, unstarch)
    #hotspot	4	https://github.com/rthurman/hotspot/archive/4.0.0.tar.gz              INSTALL
    #macs	2	macs2 2.0.10.20120913                                                     INSTALL
    #bedtools	2.17.0	http://github.com/arq5x/bedtools/archive/v2.17.0.tar.gz           COPIED TO resources/usr/bin (bedtools, bedToBam, intersectBed)
    #GCAP unversioned   https://github.com/qinqian/GCAP/                                  INSTALL
    #hotspot.py 3   https://github.com/qinqian/GCAP/gcap/glue                             COPIED to resources/usr/bin (slightly modified)

    echo "Value of genome:      '$genome'"
    echo "Value of gender:      '$gender'"
    echo "Value of ref_genome:  '$genome_ref'"
    echo "Value of chrom_sizes: '$chrom_sizes'"

    echo "* Download files..."
    #genome_ref_fn=`dx describe "$genome_ref" --name`
    dx download "$genome_ref" -o ${gentoken}.fa.gz
    gunzip ${gentoken}.fa.gz
    echo "* Genome reference fa: '${gentoken}.fa'"

    #chrom_size_fn=`dx describe "$chrom_sizes" --name`
    dx download "$chrom_sizes" -o ${gentoken}_chrom.sizes
    cat ${gentoken}_chrom.sizes | awk '{printf "%s\t0\t%s\n",$1,$2}' > ${gentoken}_chrom.bed
    echo "* chromInfo sizes and bed: '${gentoken}_chrom.sizes' '${gentoken}_chrom.bed'"

    echo "* === ${gentoken}_chrom.bed ===" 2>&1 | tee -a install.log
    cat ${gentoken}_chrom.bed | tee -a install.log
    echo "* =============================" 2>&1 | tee -a install.log

    ls -l | tee -a install.log
    ################ Installed.  Now what?
    
    #echo "* Calling peaks..."
    ##idr/bin/idr -a ${peaks_a_fn}.bed -b ${peaks_b_fn}.bed --rank-method signal.value \
    ##                     --peak-merge-method sum idr 0.20 --quiet --out-file-type bed
    #idr/bin/idr --input-file-type bed --samples peak_a.bed peaks_b_.bed --quiet
    #mv idrValues.txt ${idr_root}.bed
    #
    #echo "* Converting narrowPeak bed to bigBed..."
    #grep "^chr" ${idr_root}.bed | sort -k1,1 -k2,2n > idr_polished.bed
    #cp idr_polished.bed ${idr_root}.bb
    ##bedToBigBed idr_polished.bed -as=/usr/bin/bed6.as chromSizes.txt ${idr_root}.bb

    # Gather metrics
    meta=`echo \"DNase me\": { `
    ##          Initial parameter values: [0.10 1.00 0.20 0.50]
    #var=`grep "Initial parameter values" idr_summary.txt | awk '{printf "%s, %s, %s, %s",$4,$5,$6,$7}'`
    #var=`echo \"Initial parameter values\": $var`
    var=`echo \"About\": "nothing"`
    meta=`echo $meta $var`
    ##          Final parameter values: [0.09 0.20 0.10 0.99]
    #var=`grep "Final parameter values" idr_summary.txt | awk '{printf "%s, %s, %s, %s",$4,$5,$6,$7}'`
    #var=`echo \"Final parameter values\": $var`
    var=`echo \"Why watch\": "because it's on"`
    meta=`echo $meta, $var }`

    echo "* Upload results..."
    chrom_bed=$(dx upload ${gentoken}_chrom.bed --brief)
    log_file=$(dx upload install.log --brief)

    dx-jobutil-add-output chrom_bed "$chrom_bed" --class=file
    dx-jobutil-add-output log_file "$log_file" --class=file
    dx-jobutil-add-output metadata "$meta" --class=string

    echo "* Finished."
}
