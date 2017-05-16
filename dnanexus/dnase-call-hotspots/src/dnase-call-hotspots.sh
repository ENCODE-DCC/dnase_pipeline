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
    echo "* Value of hotspot_mappable:   '$hotspot_mappable'"

    echo "* Download files..."
    bam_root=`dx describe "$bam_to_call" --name`
    bam_root=${bam_root%.bam}
    dx download "$bam_to_call" -o ${bam_root}.bam
    echo "* bam file: '${bam_root}.bam'"

    dx download "$chrom_sizes" -o chrom.sizes
    mappable_archive=`dx describe "$hotspot_mappable" --name`
    dx download "$hotspot_mappable"
    
    hotspot_root="${bam_root}_hotspots"  # Put hotspot results into ${hotspot_root}.bed.gz, ${hotspot_root}.bb, ${hotspot_root}_count.txt, and ${hotspot_root}_SPOT.txt
    peaks_root="${bam_root}_peaks"       # Put peak results into ${peaks_root}.bed.gz, ${peaks_root}.bb, and ${peaks_root}_count.txt
    density_root="${bam_root}_density"   # Put density results into ${density_root}.bw
    allcalls_root="${bam_root}_all_calls" # Put all calls results into ${allcalls_root}.bed.gz

    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_hotspot.sh ${bam_root}.bam chrom.sizes $mappable_archive $hotspot_root $peaks_root $density_root $allcalls_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    
    echo "* Prepare metadata..."
    qc_hotspot=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        spot_score=`qc_metrics.py -n singleton -f ${hotspot_root}_SPOT.txt -k "SPOT score" --keypair "SPOT score"`
        hotspot_count=`qc_metrics.py -n singleton -f ${hotspot_root}_count.txt -k "hotspot count" --keypair "hotspot count"`
        peaks_count=`qc_metrics.py -n singleton -f ${peaks_root}_count.txt -k "peaks count" --keypair "peaks count"`
        qc_peaks=`echo $spot_score, $hotspot_count, $peaks_count`
        hotspot_count=`cat ${hotspot_root}_count.txt`
        peaks_count=`cat ${peaks_root}_count.txt`
        allcalls_count=`cat ${allcalls_root}_count.txt`
        qc_hotspot=`echo \"hotspot\": { $qc_peaks }`
    fi
    
    #### All qc to one file:
    echo "===== SPOT score ====="       > ${hotspot_root}_qc.txt
    cat ${hotspot_root}_SPOT.txt       >> ${hotspot_root}_qc.txt
    echo " "                           >> ${hotspot_root}_qc.txt
    echo "===== hotspot count ====="   >> ${hotspot_root}_qc.txt
    cat ${hotspot_root}_count.txt      >> ${hotspot_root}_qc.txt
    echo " "                           >> ${hotspot_root}_qc.txt
    echo "===== peaks count ====="     >> ${hotspot_root}_qc.txt
    cat ${peaks_root}_count.txt        >> ${hotspot_root}_qc.txt
    echo " "                           >> ${hotspot_root}_qc.txt
    echo "===== allcalls count ====="  >> ${hotspot_root}_qc.txt
    cat ${allcalls_root}_count.txt     >> ${hotspot_root}_qc.txt
    
    echo "* Upload results..."
    bed_hotspots=$(dx upload ${hotspot_root}.bed.gz   --details "{ $qc_hotspot }" --property SW="$versions" --property hotspot_count="$hotspot_count" --brief)
    bb_hotspots=$(dx upload ${hotspot_root}.bb        --details "{ $qc_hotspot }" --property SW="$versions" --property hotspot_count="$hotspot_count" --brief)
    bed_peaks=$(dx upload ${peaks_root}.bed.gz        --details "{ $qc_hotspot }" --property SW="$versions" --property peaks_count="$peaks_count" --brief)
    bb_peaks=$(dx upload ${peaks_root}.bb             --details "{ $qc_hotspot }" --property SW="$versions" --property peaks_count="$peaks_count" --brief)
    bed_allcalls=$(dx upload ${allcalls_root}.bed.gz  --details "{ $qc_hotspot }" --property SW="$versions" --property allcalls_count="$allcalls_count" --brief)
    #bw_density=$(dx upload ${density_root}.bw         --details "{ $qc_hotspot }" --property SW="$versions" --brief)
    starch_density=$(dx upload ${density_root}.starch --details "{ $qc_hotspot }" --property SW="$versions" --brief)
    hotspots_qc=$(dx upload ${hotspot_root}_qc.txt    --details "{ $qc_hotspot }" --property SW="$versions" --brief)

    dx-jobutil-add-output bed_hotspots "$bed_hotspots" --class=file
    dx-jobutil-add-output bb_hotspots "$bb_hotspots" --class=file
    dx-jobutil-add-output bed_peaks "$bed_peaks" --class=file
    dx-jobutil-add-output bb_peaks "$bb_peaks" --class=file
    dx-jobutil-add-output bed_allcalls "$bed_allcalls" --class=file
    #dx-jobutil-add-output bw_density "$bw_density" --class=file  # NOTE: this is not the density file desired.
    dx-jobutil-add-output starch_density "$starch_density" --class=file
    dx-jobutil-add-output hotspots_qc "$hotspots_qc" --class=file
    dx-jobutil-add-output metadata "{ $qc_hotspot }" --class=string
    
    echo "* Finished."
}
