#!/bin/bash
# dnase-call-hotspots.sh - Call peaks with 'hotspot' for the ENCODE DNase-seq pipeline.

main() {
    # installed pigz and gawk in dxapp.json
    alias awk='gawk'
    # hotspot and beops executables in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "* Value of bam_to_call: '$bam_to_call'"
    echo "* Value of chrom_sizes: '$chrom_sizes'"
    echo "* Value of blacklist:   '$blacklist'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_to_call" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_to_call" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    dx download "$chrom_sizes" -o chrom.sizes
    
    if [ "$blacklist" != "" ]; then
        blacklist_root=`dx describe "$blacklist" --name`
        blacklist_root=${blacklist_root%.bed.gz}
        dx download "$blacklist" -o ${blacklist}.bed.gz
        gunzip ${blacklist}.bed.gz
        blacklist_file="${blacklist}.bed"
        echo "* blacklist file: '$blacklist_file'"
    else
        blacklist_file=""
    fi
    
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_hotspot.sh ${bam_root}.bam chrom.sizes $blacklist_file
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    hotspot_root="${bam_root}_hotspots"
    peaks_root="${bam_root}_peaks"
    density_root="${bam_root}_density"
    ### Temporary for debugging
    cutcounts_root="${bam_root}_cutcounts"
    allcalls_root="${bam_root}_allcalls"
    ### Temporary for debugging
    
    echo "* Compressing bed file..."
    set -x
    pigz ${hotspot_root}.bed
    pigz ${peaks_root}.bed
    set +x

    echo "* Prepare metadata..."
    qc_hotspot=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        ###    qc_hotspot=`qc_metrics.py -n hotspot -f ${bam_root}_hotspot_out.txt`
        hotspot_count=`qc_metrics.py -n singleton -f ${hotspot_root}_count.txt -k "hotspot count" --keypair "hotspot count"`
        peaks_count=`qc_metrics.py -n singleton -f ${peaks_root}_count.txt -k "peaks count" --keypair "peaks count"`
        qc_peaks=`echo $hotspot_count, $peaks_count`
        hotspot_count=`cat ${hotspot_root}_count.txt`
        peaks_count=`cat ${peaks_root}_count.txt`
        ### Temporary for debugging
        cut_count=`qc_metrics.py -n singleton -f ${cutcounts_root}_count.txt -k "cut count" --keypair "cut count"`  
        allcalls_count=`qc_metrics.py -n singleton -f ${allcalls_root}_count.txt -k "all calls" --keypair "all calls"`  
        qc_peaks=`echo $qc_peaks, $cut_count, $allcalls_count`
        cut_count=`cat ${cutcounts_root}_count.txt`  
        allcalls_count=`cat ${allcalls_root}_count.txt`  
        ### Temporary for debugging
        #qc_hotspot=`echo $qc_hotspot, \"peak_counts\": { $qc_spots, $qc_regions }`
        qc_hotspot=`echo \"peak_counts\": { $qc_peaks }`
    fi
    #### All qc to one file:
    echo "===== hotspot count ====="    > ${hotspot_root}_qc.txt
    cat ${hotspot_root}_count.txt      >> ${hotspot_root}_qc.txt
    echo " "                           >> ${hotspot_root}_qc.txt
    echo "===== peaks count ====="     >> ${hotspot_root}_qc.txt
    cat ${peaks_root}_count.txt        >> ${hotspot_root}_qc.txt
    echo " "                           >> ${hotspot_root}_qc.txt    ### Temporary for debugging
    echo "===== allcalls count ====="  >> ${hotspot_root}_qc.txt    ### Temporary for debugging
    cat ${allcalls_root}_count.txt     >> ${hotspot_root}_qc.txt    ### Temporary for debugging
    echo " "                           >> ${hotspot_root}_qc.txt    ### Temporary for debugging
    echo "===== cut count ====="       >> ${hotspot_root}_qc.txt    ### Temporary for debugging
    cat ${cutcounts_root}_count.txt    >> ${hotspot_root}_qc.txt    ### Temporary for debugging
    
    echo "* Upload results..."
    bed_hotspots=$(dx upload ${hotspot_root}.bed.gz  --details "{ $qc_hotspot }" --property SW="$versions" --property hotspot_count="$hotspot_count" --brief)
    bb_hotspots=$(dx upload ${hotspot_root}.bb       --details "{ $qc_hotspot }" --property SW="$versions" --property hotspot_count="$hotspot_count" --brief)
    bed_peaks=$(dx upload ${peaks_root}.bed.gz       --details "{ $qc_hotspot }" --property SW="$versions" --property peaks_count="$peaks_count" --brief)
    bb_peaks=$(dx upload ${peaks_root}.bb            --details "{ $qc_hotspot }" --property SW="$versions" --property peaks_count="$peaks_count" --brief)
    bw_density=$(dx upload ${density_root}.bw        --details "{ $qc_hotspot }" --property SW="$versions" --brief)
    hotspots_qc=$(dx upload ${hotspot_root}_qc.txt   --details "{ $qc_hotspot }" --property SW="$versions" --brief)
      ### Temporary for debugging
    starch_cutcounts=$(dx upload ${cutcounts_root}.starch --details "{ $qc_hotspot }" --property SW="$versions" --property cut_count="$cut_count" --brief)
      ### Temporary for debugging
    starch_allcalls=$(dx upload ${allcalls_root}.starch   --details "{ $qc_hotspot }" --property SW="$versions" --property allcalls_count="$allcalls_count" --brief)

    dx-jobutil-add-output bed_hotspots "$bed_hotspots" --class=file
    dx-jobutil-add-output bb_hotspots "$bb_hotspots" --class=file
    dx-jobutil-add-output bed_peaks "$bed_peaks" --class=file
    dx-jobutil-add-output bb_peaks "$bb_peaks" --class=file
    dx-jobutil-add-output bw_density "$bw_density" --class=file
    echo "* Add hotspots_qc to output..."
    dx-jobutil-add-output hotspots_qc "$hotspots_qc" --class=file
    echo "* Add starch_cutcounts to output..."
    dx-jobutil-add-output starch_cutcounts "$starch_cutcounts" --class=file  ### Temporary for debugging
    echo "* Add starch_allcalls to output..."
    dx-jobutil-add-output starch_allcalls "$starch_allcalls" --class=file  ### Temporary for debugging
    echo "* Add metadata to output..."
    dx-jobutil-add-output metadata "{ $qc_hotspot }" --class=string

    echo "* Finished."
}
