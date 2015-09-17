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
- dnase-align-bwe-pe  - Takes two (paired-end) gzipped fastq file sets and a bwa built index of a gender specific
                        reference genome (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns
                        with bwa and returns a single bam file containing all reads mapped or otherwise and a QC metrics file.  
                        Used in ENCODE to align technical replicates.
- dnase-align-bwe-se  - Takes one (single-end) gzipped fastq file set and a bwa built index of a gender specific 
                        reference genome (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns
                        with bwa and returns a single bam file containing all reads mapped or otherwise and a QC metrics file.
                        Used in ENCODE to align technical replicates.
- dnase-filter-pe     - Takes one or more paired-end bams produced by '`dnase-align-bwe-pe`' and merges then filters to 
                        include only high quality mapped reads.  This step produces a filtered bam and a QC metrics file.
- dnase-filter-se     - Takes one or more bams produced by '`dnase-align-bwe-se`' (or '`dnase-align-bwe-pe`') 
                        and merges then filters to include only high quality mapped reads.  This step produces a 
                        filtered bam and a QC metrics file.
- dnase-eval-bam-pe   - Takes one paired-end bam produced by '`dnase-filter-pe`', excludes mappings to chrM, sub-samples
                        reads and generates evaluation data. This step produces a sample bam and a QC metrics file.
- dnase-eval-bam-se   - Takes one single-ended bam produced by '`dnase-filter-se`', excludes mappings to chrM, sub-samples
                        reads and generates evaluation data. This step produces a sample bam and a QC metrics file.
- dnase-hotspot-qc    - Takes the bam filtered from either paired or single ended datasets.
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
INPUTS:  read1.fq.gz            n*tech-pe.bam(b)  filtered-pe.bam(c)  filtered-pe.bam(c)  filtered-pe.bam(c)   2*filtered-pe.bam(c)
         read2.fq.gz                  |                   |                  |            chrom.sizes          2*hs_narrowPeak.bb(d)
         bwa_index.tgz(a)             |                   |                  |                   |             2*hs_signal.bw(e)
              |                       |                   |                  |                   |                      |
              V                       V                   V                  V                   V                      V
STEPS:   dnase-align-bwa-pe ==> bam-filter-pe ==> dnase-eval-bam ==> dnase-hotspot-qc ==> dnase-call-hotspots ==> dnase-pool-bioreps
              |                       |                   |                  |                   |                      |
              V                       V                   V                  V                   V                      V
OUTPUTS: tech-pe.bam(b)         filtered-pe.bam(c)  sampled-pe.bam     sample_5M.bam      hs_narrowPeak.bb,bed(d)  pooled.bam
                                filtered-pe-qc.txt  sampled-pe-qc.txt  sample_5M_qc.txt   hs_broadPeak.bb,bed      merged.bb
                                                                                          hotspot_signal.bw(e)     merged.bed
                                                                                          hotspot_qc.txt           pooled_qc.txt
```

---------
*Single-end pipeline:*
```
INPUTS:  reads.fq.gz            n*tech-se.bam(f)  filtered-se.bam(g)  filtered-se.bam(g)  filtered-se.bam(g)      2*filtered-pe.bam(g)
         bwa_index.tgz(a)       or tech-pe.bam(b)         |                  |            chrom.sizes             2*hs_narrowPeak.bb(h)
              |                       |                   |                  |                   |                2*hs_signal.bw(i)
              |                       |                   |                  |                   |                      |
              V                       V                   V                  V                   V                      V
STEPS:   dnase-align-bwa-se ==> bam-filter-se ==> dnase-eval-bam ==> dnase-hotspot-qc ==> dnase-call-hotspots ==> dnase-pool-bioreps
              |                       |                   |                  |                   |                      |
              V                       V                   V                  V                   V                      V
OUTPUTS: tech-se.bam(f)         filtered-se.bam(g)  sampled-se.bam     sample_5M.bam      hs_narrowPeak.bb,bed(h)  pooled.bam(q)
                                filtered-se-qc.txt  sampled-se-qc.txt  no_chrM-se-qc.txt  hs_broadPeak.bb,bed      merged.bb
                                                                                          hotspot_signal.bw(i)     merged.bed
                                                                                          hotspot_qc.txt           pooled_qc.txt
```

