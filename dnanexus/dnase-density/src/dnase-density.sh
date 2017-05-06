#!/bin/bash
# dnase-eval-bam.sh - Evaluates sample of bam for the ENCODE DNase-seq pipeline.

main() {
    # executables in resources/usr/bin
    set +x
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi
 
    echo "* Value of bam_filtered:  '$bam_filtered'"
    echo "* Value of chrom_sizes:   '$chrom_sizes'"

    echo "* Download files..."
    # expecting *_bwa_biorep_filtered.bam
    bam_input_root=`dx describe "$bam_filtered" --name`
    # Better to leave the whole suffix!
    bam_input_root=${bam_input_root%.bam}
    dx download "$bam_filtered" -o ${bam_input_root}.bam
    echo "* bam file: '${bam_input_root}.bam'"

    dx download "$chrom_sizes" -o chrom.sizes
    
    density_root="${bam_input_root}_normalized_density"
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    dnase_density.sh ${bam_input_root}.bam chrom.sizes $density_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="

    echo "* Prepare metadata for density files..."
    qc_density=''
    # TODO: Any qc?
    #if [ -f /usr/bin/qc_metrics.py ]; then
    #    qc_density=`qc_metrics.py -n edwBamStats -f ${bam_sample_root}_edwBamStats.txt`
    #    meta=`qc_metrics.py -n phantompeaktools_spp -f ${bam_sample_root}_spp_out.txt`
    #    qc_density=`echo $qc_density, $meta`
    #fi
    ## All qc to one file per target file:
    #echo "===== edwBamStats ====="           > ${density_root}_qc.txt
    #cat ${bam_sample_root}_edwBamStats.txt  >> ${density_root}_qc.txt
    #echo " "                                >> ${density_root}_qc.txt
    #echo "===== phantompeaktools spp =====" >> ${density_root}_qc.txt
    #cat ${bam_sample_root}_spp_out.txt      >> ${density_root}_qc.txt
        
    echo "* Upload results..."
    norm_density_bw=$(dx upload ${density_root}.bw --details "{ $qc_density }" --property SW="$versions" --brief)
    norm_density_starch=$(dx upload ${density_root}.starch  --details "{ $qc_density }" --property SW="$versions" --brief)
    #norm_density_qc=$(dx upload ${density_root}_qc.txt --details "{ $qc_density }" --property SW="$versions" --brief)

    dx-jobutil-add-output norm_density_bw "$norm_density_bw" --class=file
    dx-jobutil-add-output norm_density_starch "$norm_density_starch" --class=file
    #dx-jobutil-add-output norm_density_qc "$norm_density_qc" --class=file

    dx-jobutil-add-output metadata "{ $qc_density }" --class=string

    echo "* Finished."
}
