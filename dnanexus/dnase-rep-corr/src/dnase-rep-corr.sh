#!/bin/bash
# dnase-idr.sh

main() {
    # installed r-base-core in dxapp.json
    # chromCor.Rscript, bigWigToWig and beops starch executables in resources/usr/bin
    
    # If available, will print tool versions to stderr and json string to stdout
    versions=''
    if [ -f /usr/bin/tool_versions.py ]; then 
        versions=`tool_versions.py --dxjson dnanexus-executable.json`
    fi

    echo "Value of density_a:   '$density_a'"
    echo "Value of density_b:   '$density_b'"

    echo "* Download files..."
    format="unknown"
    density_a_file=`dx describe "$density_a" --name`
    density_a_root=${density_a_file%.bw}
    density_a_root=${density_a_root%.bigWig}
    if [[ "$density_a_root" != "$density_a_file" ]]; then
        format="bigWig"
    else
        density_a_root=${density_a_file%.starch}
        if [[ "$density_a_root" != "$density_a_file" ]]; then
            format="starch"
        fi
    fi
    if [ "$format" == "unknown" ]; then
        echo "* ERROR: could not detect format of '$density_a_file'.  Only 'starch' and 'bigWig' are supported."
        exit 1
    fi
    # ex: ENCSR691MQJ_rep2_1_se_bwa_biorep_filtered_normalized_density.bw => ENCSR691MQJ_rep2_1
    density_a_root=${density_a_root%_density}
    density_a_root=${density_a_root%_normalized}
    density_a_root=${density_a_root%_filtered}
    density_a_root=${density_a_root%_biorep}
    density_a_root=${density_a_root%_bwa}
    density_a_root=${density_a_root%_se}
    density_a_root=${density_a_root%_pe}
    density_a_root=${density_a_root%_spe}
    if [ "$format" == "bigWig" ]; then
        echo "* Download density file: '${density_a_root}.bw'"
        dx download "$density_a" -o ${density_a_root}.bw
        echo "* Convert 'bigWig' to 'starch' file: '${density_a_root}.starch'"
        set -x
        bigWigToWig ${density_a_root}.bw stdout | grep -v \#bedGraph | awk -v OFS='\t' '{print $1,$2,$3,".",$4}' | \
                                                                                    starch - > ${density_a_root}.starch
        set +x
    else
        echo "* Download density file: '${density_a_root}.starch'"
        dx download "$density_a" -o ${density_a_root}.starch
    fi

    format="unknown"
    density_b_file=`dx describe "$density_b" --name`
    density_b_root=${density_b_file%.bw}
    density_b_root=${density_b_root%.bigWig}
    if [[ "$density_b_root" != "$density_b_file" ]]; then
        format="bigWig"
    else
        density_b_root=${density_b_file%.starch}
        if [[ "$density_b_root" != "$density_b_file" ]]; then
            format="starch"
        fi
    fi
    if [ "$format" == "unknown" ]; then
        echo "* ERROR: could not detect format of '$density_b_file'.  Only 'starch' and 'bigWig' are supported."
        exit 1
    fi
    density_b_root=${density_b_root%_density}
    density_b_root=${density_b_root%_normalized}
    density_b_root=${density_b_root%_filtered}
    density_b_root=${density_b_root%_biorep}
    density_b_root=${density_b_root%_bwa}
    density_b_root=${density_b_root%_se}
    density_b_root=${density_b_root%_pe}
    density_b_root=${density_b_root%_spe}
    if [ "$format" == "bigWig" ]; then
        echo "* Download density file: '${density_b_root}.bw'"
        dx download "$density_b" -o ${density_b_root}.bw
        echo "* Convert 'bigWig' to 'starch' file: '${density_b_root}.starch'"
        set -x
        bigWigToWig ${density_b_root}.bw stdout | grep -v \#bedGraph | awk -v OFS='\t' '{print $1,$2,$3,".",$4}' | \
                                                                                    starch - > ${density_b_root}.starch
        set +x
    else
        echo "* Download density file: '${density_b_root}.starch'"
        dx download "$density_b" -o ${density_b_root}.starch
    fi

    corr_root="${density_a_root}_to_${density_b_root}_normalized_density_corr"
    echo "* Correlation root: '"$corr_root"'"

    # DX/ENCODE independent script is found in resources/usr/bin
    echo "* ===== Calling DNAnexus and ENCODE independent script... ====="
    set -x
    cp /usr/bin/chromCor.Rscript .
    dnase_rep_corr.sh ${density_a_root}.starch ${density_b_root}.starch 5 $corr_root
    set +x
    echo "* ===== Returned from dnanexus and encodeD independent script ====="
    
    echo "* Prepare metadata..."
    qc_stats=''
    if [ -f /usr/bin/qc_metrics.py ]; then
        qc_stats=`qc_metrics.py -n rep_corr -f ${corr_root}.txt`
    fi
    
    echo "* Upload results..."
    corr_txt=$(dx upload ${corr_root}.txt --details="{ $qc_stats }" --property SW="$versions" --brief)

    dx-jobutil-add-output corr_txt "$corr_txt" --class=file
    dx-jobutil-add-output metadata "{ $qc_stats }" --class=string

    echo "* Finished."
}
