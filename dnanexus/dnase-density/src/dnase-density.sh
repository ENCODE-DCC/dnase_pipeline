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
    density_count=""
    # qc anyone?  Not much value.
    if [ -f ${density_root}_count.txt ]; then
        density_count=`cat ${density_root}_count.txt`
        ## All qc to one file per target file:
        #echo "===== density intervals ====="  > ${density_root}_qc.txt
        #cat ${density_root}_count.txt        >> ${density_root}_qc.txt
    fi
        
    echo "* Upload results..."
    normalized_bw=$(dx upload ${density_root}.bw --property intervals=$density_count --property SW="$versions" --brief)
    normalized_starch=$(dx upload ${density_root}.starch  --property intervals=$density_count --property SW="$versions" --brief)
    #normalized_qc=$(dx upload ${density_root}_qc.txt --property intervals=$density_count --property SW="$versions" --brief)

    dx-jobutil-add-output normalized_bw "$normalized_bw" --class=file
    dx-jobutil-add-output normalized_starch "$normalized_starch" --class=file
    #dx-jobutil-add-output norm_density_qc "$norm_density_qc" --class=file

    echo "* Finished."
}
