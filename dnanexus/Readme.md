# ENCODE dnase-seq-pipeline (DNAnexus Platform App)

This folder contains the dnanexis applets used in the ENCODE dnase-seq-pipeline. Most applets run 
bash scripts with straight-forward command-lines for third party tools.  The command-lines can readily 
be repurposed for a non-dnanexus environment.

## Applets (aka steps):
*Preparation:* Typically this step is only one once per gender specific reference genome (e.g. GCRh37/hg19 female).  
               The result of the preparation step is used in all subsequent pipeline runs for the same genome.
- index-bwa - Takes a gender specific genome reference gzipped fasta file (e.g. GCRh37/hg19 female) and produces a 
              bwa index in the form of a gzipped tar file.

*Normal pipline:* The normal pipeline is actually 2 separate pipelines, as paired-end and single end fastqs are handled
                differently.  The '-pe' and '-se' alignment steps encapsulate the differences. 
- align-bwe-pe     - Takes two (paired-end) gzipped fastqs and a bwa build index of a gender specific reference genome
                     (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns with bwa
                     and returns a single bam file.
- align-bwe-se     - Takes one (single-end) gzipped fastq and a bwa build index of a gender specific reference genome
                     (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns with bwa
                     and returns a single bam file containing all reads mapped or otherwise.
- bam-filter-pe    - Takes the bwa aligned bam produced by '`align-bwa-pe`' and filters to include only high quality
                     mapped reads.  This step produces a filtered bam and an one that further excludes chrM mappings,
                     as well as two files of QC metrics.
- bam-filter-se    - Takes the bwa aligned bam produced by '`align-bwa-se`' and filters to include only high quality
                     mapped reads.  This step produces a filtered bam and an one that further excludes chrM mappings,
                     as well as two files of QC metrics.
- sample-hotspots  - Takes the bam filtered and without chrM mappings from either paired or single ended datasets.
                     This step runs hotspot on a 5 million read sample of the bam in order to produce QC metrics.
- call-hotspots    - Taks a filtered bam and runs hotspots.  This step produces narrowPeaks (hotspots) and broadPeaks
                     (regions) in both bed and bigBed format.  The step also produces signals in bigWig format and
                     a hotspot QC metrics file.  Note that '`call-hotspots`' is run on both the filtered non-chrM bams 
                     of individual replicates, and on the bam produced by merging two replicates.
- merge-replicates - Takes filtered bams, narrowPeaks (hotspots) and signals from two replicates and produces a merge 
                     bam, merged hotspots and two QC metrics files comparing the two replicates.

---------
## Flow
*Prepartion step:*
```
INPUTS:  genome-gender.fa.gz
             |
             V
STEPS:   index-bwa
             |
             V
OUTPUTS: bwa_index.tgz(a)
```

---------
*Paired-end pipeline:*
```
INPUTS:  read1.fq.gz       bwa-pe.bam(b)      no_chrM-pe.bam(d)    no_chrM-pe.bam(d)   2*filtered-pe.bam(c)   merged.bam(g)
         read2.fq.gz            |                   |              chrom.sizes         2*hs_narrowPeak.bb(e)  chrom.sizes 
         bwa_index.tgz(a)       |                   |                    |             2*hs_signal.bw(f)          |
              |                 |                   |                    |                  |                     |
              V                 V                   V                    V                  V                     V
STEPS:   align-bwa-pe ===> bam-filter-pe ===> sample_hotspots ===> call-hotspots ===>  merge-replicates ===> call-hotspots
              |                 |                   |                    |                  |                     |
              V                 V                   V                    V                  V                     V
OUTPUTS: bwa-pe.bam(b)     filtered-pe.bam(c) hotspot_5M_qc.txt hs_narrowPeak.bb,bed(e) merged.bam(g)        pooled_hs_narrowPeak.bb,bed
                           no_chrM-pe.bam(d)  no_chrM_5M.bam    hs_broadPeak.bb,bed     merged.bb            pooled_hs_broadPeak.bb,bed
                           filtered-qc.txt                      hotspot_signal.bw(f)    merged.bed           pooled_hotspot_signal.bw
                           filtered-qc-full.txt                 hotspot_qc.txt          signal_corr_qc.txt   pooled_hotspot_qc.txt
                                                                                        peaks_overlap_qc.txt
```

---------
*Single-end pipeline:*
```
INPUTS:  reads.fq.gz       bwa-se.bam(h)      no_chrM-pe.bam(j)    no_chrM-pe.bam(j)   2*filtered-pe.bam(i)     merged.bam(m)
         bwa_index.tgz(a)       |                   |              chrom.sizes         2*hs_narrowPeak.bb(k)    chrom.sizes 
              |                 |                   |                    |             2*hs_signal.bw(l)          |
              |                 |                   |                    |                  |                     |
              V                 V                   V                    V                  V                     V
STEPS:   align-bwa-se ===> bam-filter-se ===> sample_hotspots ===> call-hotspots ===>  merge-replicates ===> call-hotspots
              |                 |                   |                    |                  |                     |
              V                 V                   V                    V                  V                     V
OUTPUTS: bwa-se.bam(h)     filtered-se.bam(i) hotspot_5M_qc.txt hs_narrowPeak.bb,bed(k) merged.bam(m)        pooled_hs_narrowPeak.bb,bed
                           no_chrM-se.bam(j)  no_chrM_5M.bam    hs_broadPeak.bb,bed     merged.bb            pooled_hs_broadPeak.bb,bed
                           filtered-qc.txt                      hotspot_signal.bw(l)    merged.bed           pooled_hotspot_signal.bw
                           filtered-qc-full.txt                 hotspot_qc.txt          signal_corr_qc.txt   pooled_hotspot_qc.txt
                                                                                        peaks_overlap_qc.txt
```


                                                         
