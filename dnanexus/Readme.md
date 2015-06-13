# ENCODE dnase-seq-pipeline (DNAnexus Platform App)

This folder contains the dnanexis applets used in the ENCODE dnase-seq-pipeline. Most applets run 
bash scripts with straight-forward command-lines for third party tools.  The command-lines can readily 
be repurposed for a non-dnanexus environment.

NEW PLAN:
1) align each tech_rep
2) Combine tech_reps AND fileter.  Set output param filtered bam read count
   - save filtered bam
3) Subsample to size with top limit and sister bio_rep filtered bam read count (maintains derived_from)
   - save sized bam (can it rename large bam?)
4) QA and no_chrM polishing
   - save no_chrM bam
5) call_hotspots
6) sample_hotspots
7) merge_replicates
8) call_hotspots again.

## Applets (aka steps):
*Preparation:* Typically this step is only one once per gender specific reference genome (e.g. GCRh37/hg19 female).  
               The result of the preparation step is used in all subsequent pipeline runs for the same genome.
- index-bwa - Takes a gender specific genome reference gzipped fasta file (e.g. GCRh37/hg19 female) and produces a 
              bwa index in the form of a gzipped tar file.

*Normal pipline:* The normal pipeline is actually 2 separate pipelines, as paired-end and single end fastqs are handled
                differently.  The '-pe' and '-se' alignment steps encapsulate the differences. 
- dnase-align-bwe-pe  - Takes two (paired-end) gzipped fastqs and a bwa build index of a gender specific reference genome
                        (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns with bwa
                        and returns a single bam file.  Used in ENCODE to align technical replicates.
- dnase-align-bwe-se  - Takes one (single-end) gzipped fastq and a bwa build index of a gender specific reference genome
                        (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns with bwa
                        and returns a single bam file containing all reads mapped or otherwise.
                        Used in ENCODE to align technical replicates.
- dnase-merge-bams    - Takes two or more bam produced by '`dnase-align-bwe-se`' or '`dnase-align-bwe-pe`' and returns
                        a single bam.  Used in ENCODE to combine multiple technical replicates for a single biological
                        repliecate.
- dnase-filter-pe     - Takes one paired-end bam produced by '`dnase-merge-bams`' and filters to include only high 
                        quality mapped reads.  This step produces a filtered bam and two files of QC metrics.
- dnase-filter-se     - Takes one single-ended bam produced by '`dnase-merge-bams`' and filters to include only high 
                        quality mapped reads.  This step produces a filtered bam and two files of QC metrics.
- dnase-size-bam      - Takes a bam produced by '`dnase-filter-pe`' or '`dnase-filter-se`' and returns
                        a bam sized to a target number of reads.  Used in ENCODE to ensure two bams from biological 
                        replicates will produce comparable results when hotspots are called.
- dnase-eval-bam-pe   - Takes one paired-end bam produced by '`dnase-size-bam`', excludes mappings to chrM, sub-samples
                        reads and generates evaluation data. This step produces the bam without chrM, a smaller sample
                        bam and four files of QC metrics.
- dnase-eval-bam-se   - Takes one single-ended bam produced by '`dnase-size-bam`', excludes mappings to chrM, sub-samples
                        reads and generates evaluation data. This step produces the bam without chrM, a smaller sample
                        bam and four files of QC metrics.
- dnase-hotspot-qc    - Takes the bam filtered and without chrM mappings from either paired or single ended datasets.
                        This step runs hotspot on a 5 million read sample of the bam in order to produce QC metrics.
- dnase-call-hotspots - Takes a filtered bam and runs hotspots.  This step produces narrowPeaks (hotspots) and broadPeaks
                        (regions) in both bed and bigBed format.  The step also produces signals in bigWig format and
                        a hotspot QC metrics file.  Note that '`call-hotspots`' is run on both the filtered non-chrM bams 
                        of individual replicates, and on the bam produced by merging two replicates.
- dnase-pool-bioreps  - Takes filtered bams, narrowPeaks (hotspots) and signals from two replicates and produces a pooled 
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
INPUTS:  read1.fq.gz            n*tech-pe.bam(b)    bio-rep-pe.bam(c)  filtered-pe.bam(d)  sized-pe.bam(e)   no_chrM-pe.bam(f)   no_chrM-pe.bam(f)      2*sized-pe.bam(e)       pooled.bam(i)
         read2.fq.gz                  |                   |                    |                 |                  |            chrom.sizes            2*hs_narrowPeak.bb(g)   chrom.sizes 
         bwa_index.tgz(a)             |                   |                    |                 |                  |                  |                2*hs_signal.bw(h)           |
              |                       |                   |                    |                 |                  |                  |                       |                    |
              V                       V                   V                    V                 V                  V                  V                       V                    V
STEPS:   dnase-align-bwa-pe ==> dnase-merge-bams ==> bam-filter-pe ==> dnase-size-bam ==> dnase-eval-bam ==> dnase-hotspot-qc ==> dnase-call-hotspots ==> dnase-pool-bioreps ==> dnase-call-hotspots
              |                       |                    |                   |                 |                  |                  |                       |                    |
              V                       V                    V                   V                 V                  V                  V                       V                    V
OUTPUTS: tech-pe.bam(b)         bio-rep-pe.bam(c)  filtered-pe.bam(d)    sized-pe.bam(e)  no_chrM-pe.bam(f)  sample_5M.bam      hs_narrowPeak.bb,bed(g)   pooled.bam(i)   pooled_hs_narrowPeak.bb,bed
                                bio-rep-pe-qc.txt  filtered-pe-qc.txt    sized-pe-qc.txt  no_chrM-pe-qc.txt  sample_5M_qc.txt   hs_broadPeak.bb,bed       merged.bb       pooled_hs_broadPeak.bb,bed
                                                                                          sampled-pe.bam                        hotspot_signal.bw(h)      merged.bed      pooled_hotspot_signal.bw
                                                                                          sampled-pe-qc.txt                     hotspot_qc.txt            pooled_qc.txt   pooled_hotspot_qc.txt
```

---------
*Single-end pipeline:*
```
INPUTS:  read1.fq.gz            n*tech-se.bam(j)    bio-rep-se.bam(k)  filtered-se.bam(l)  sized-se.bam(m)   no_chrM-se.bam(n)   no_chrM-pe.bam(n)      2*sized-se.bam(m)       pooled.bam(q)
         read2.fq.gz            or tech-pe.bam(a)         |                    |                 |                  |            chrom.sizes            2*hs_narrowPeak.bb(o)   chrom.sizes 
         bwa_index.tgz(a)             |                   |                    |                 |                  |                  |                2*hs_signal.bw(p)           |
              |                       |                   |                    |                 |                  |                  |                       |                    |
              V                       V                   V                    V                 V                  V                  V                       V                    V
STEPS:   dnase-align-bwa-se ==> dnase-merge-bams ==> bam-filter-se ==> dnase-size-bam ==> dnase-eval-bam ==> dnase-hotspot-qc ==> dnase-call-hotspots ==> dnase-pool-bioreps ===> dnase-call-hotspots
              |                       |                    |                   |                 |                  |                  |                       |                    |
              V                       V                    V                   V                 V                  V                  V                       V                    V
OUTPUTS: tech-se.bam(b)         bio-rep-se.bam(k)  filtered-se.bam(l)    sized-se.bam(m)  no_chrM-se.bam(n)  sample_5M.bam      hs_narrowPeak.bb,bed(o)   pooled.bam(q)   pooled_hs_narrowPeak.bb,bed
                                bio-rep-se-qc.txt  filtered-se-qc.txt    sized-se-qc.txt  no_chrM-se-qc.txt  sample_5M_qc.txt   hs_broadPeak.bb,bed       merged.bb       pooled_hs_broadPeak.bb,bed
                                                                                          sampled-se.bam                        hotspot_signal.bw(p)      merged.bed      pooled_hotspot_signal.bw
                                                                                          sampled-se-qc.txt                     hotspot_qc.txt            pooled_qc.txt   pooled_hotspot_qc.txt
```




## ------- Older design now replaced -------
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


                                                         
