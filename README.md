dnase-pipeline
==============

ENCODE / DNA nexus pipeline for DNase-seq experiments.

Steps:
- Alignment by bwa
- Filtering of bams and QC analysis
- Peak/hotspot calling by hotspot
- Sampled hotspot for QC analysis
- Merged replicates for QC analysis
- Peak/hotspot calls on merged replicates
