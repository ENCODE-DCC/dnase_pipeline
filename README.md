dnase-pipeline
==============

ENCODE Uniform Processing Pipeline for DNaseHS experiments.

Enables running on DNANexus but each step is available as a stand alone (DNANexus and ENCODE independent) script.

Steps:
* Preparation: bwa indexing and hotspot mappability/blacklist preparation (run once on reference files)
  * stand-alone script: dnanexus/dnase-index-bwa/resources/usr/bin/dnase_index_bwa.sh
* Alignment by bwa
  * paired-end script: dnanexus/dnase-align-bwa-pe/resources/usr/bin/dnase_align_bwa_pe.sh
  * single-end script: dnanexus/dnase-align-bwa-se/resources/usr/bin/dnase_align_bwa_se.sh
* Merging and Filtering of bams and QC analysis
  * paired-end script: dnanexus/dnase-filter-pe/resources/usr/bin/dnase_filter_pe.sh
  * single-end script: dnanexus/dnase-filter-se/resources/usr/bin/dnase_filter_se.sh
* BAM evaluation script
  * script: dnanexus/dnase-eval-bam/resources/usr/bin/dnase_eval_bam.sh
* Peak/hotspot calling by hotspot
  * script: dnanexus/dnase-call-hotspots/resources/usr/bin/dnase_hotspot.sh
* Replicate Concordance
  * script: dnanexus/dnase-rep-corr/resources/usr/bin/dnase_rep_corr.sh

Further detail and basic workflows can be seen in dnanexus/Readme.md.  Most dependencies for the above scripts 
are contained in the same bin directories as the scripts.

There are tools to help set up a DNANexus project with these scripts:
* dnanexus/build_applets will build the applets into the named project
  * This script will copy in common tools used by all steps and found in the dnanexus/tools directory
  * Note that there are separate paired/single end versions of the pipelines and common steps be created distinct versions.
* dnanexus/dnaseLaunch.py can be used to build template workflows.  However, it is dependent upon 
  * https://github.com/ENCODE-DCC/dxencode/blob/master/template.py


