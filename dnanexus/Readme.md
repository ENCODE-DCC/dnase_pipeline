# ENCODE Uniform Processing Pipeline for DNaseHS experiments (dnase-seq-pipeline)

This folder contains the dnanexus applets used in the ENCODE dnase-seq-pipeline. Most applets run 
bash scripts with straight-forward command-lines for third party tools.  The command-lines can readily 
be repurposed for a non-dnanexus environment.

## Applets (aka steps):
*Preparation:* Typically this step is only run once per reference genome (e.g. GCRh37/hg19).  
               The result of the preparation step is used in all subsequent pipeline runs for the same genome.
- dnase-index-bwa - Takes a genome reference gzipped fasta file (e.g. GCRh37/hg19), a DNase mappability file (optional) 
                    specific to a read size, a genome blacklist file (optional) and produces a bwa index in the form of 
                    a gzipped tar file, and a mappability gzipped tar file.
*Normal pipline:* The normal pipeline is actually 2 separate pipelines, as paired-end and single end fastqs are handled
                differently.  The '-pe' and '-se' alignment steps encapsulate the differences. 
- dnase-align-bwe-pe  - Takes two (paired-end) gzipped fastq file sets and a bwa built index of a reference genome 
                        (e.g. GCRh37/hg19) in the form of a gzipped tar archive.  This step aligns with bwa and returns 
                        a single bam file containing all reads mapped or otherwise and a QC metrics file.  
                        Used in ENCODE to align technical replicates.
- dnase-align-bwe-se  - Takes one (single-end) gzipped fastq file set and a bwa built index of a reference genome 
                        (e.g. GCRh37/hg19 female) in the form of a gzipped tar archive.  This step aligns with bwa and 
                        returns a single bam file containing all reads mapped or otherwise and a QC metrics file.
                        Used in ENCODE to align technical replicates.
- dnase-filter-pe     - Takes one or more paired-end bams produced by '`dnase-align-bwe-pe`' and merges then filters to 
                        include only high quality mapped reads.  This step produces a filtered bam and a QC metrics file.
                        Used in ENCODE to merge technical replicates into a single biological replicate filtered bam file.
- dnase-filter-se     - Takes one or more bams produced by '`dnase-align-bwe-se`' (or '`dnase-align-bwe-pe`') 
                        and merges then filters to include only high quality mapped reads.  This step produces a 
                        filtered bam and a QC metrics file.
                        Used in ENCODE to merge technical replicates into a single biological replicate filtered bam file.
- dnase-eval-bam      - Takes one paired-end bam produced by '`dnase-filter-pe`' or '`dnase-filter-se`', excludes mappings 
                        to chrM, sub-samples reads and generates evaluation data. This step produces a sample bam and a 
                        QC metrics file.
                        Used in ENCODE for QC metrics only.  The subsampled bam file is not saved.
- dnase-call-hotspots - Takes a filtered bam produced by '`dnase-filter-pe`' or '`dnase-filter-se`' and runs hotspots.
                        This step produces narrowPeaks (hotspots) and broadPeaks (regions) in both bed and bigBed format.  
                        The step also produces signals in bigWig format and a hotspot QC metrics file.  
                        Used in ENCODE to produce the target results of DNase experiments.
- dnase-rep-corr      - Takes two density (*.bw) files produced by two runs of '`dnase-call-hotspots`' on separate replicates
                        of the same experiment.  This step produces a single QC metrics file.
                        Used in ENCODE to ensure the reproducibility of results with seperate biological replicates.

---------
## Flow
*Prepartion step:*
```
INPUTS:  genome.fa.gz, genome-blacklist.bed.gz, mappability_only.bed.gz 
             |
             V
STEPS:   index-bwa
             |
             V
OUTPUTS: bwa_index.tgz(a), mappability.tgz(b)
```

---------
*Paired-end pipeline:*
```
INPUTS:  read1.fq.gz            n*tech-pe.bam(c)  filtered-pe.bam(d)  filtered-pe.bam(d)   2*hotspot_signal.bw(f)
         read2.fq.gz                  |                   |           mappability.tgz(b)          |
         bwa_index.tgz(a)             |                   |           chrom.sizes                 |
              |                       |                   |                  |                    |
              V                       V                   V                  V                    V
STEPS:   dnase-align-bwa-pe ==> bam-filter-pe ==> dnase-eval-bam ==> dnase-call-hotspots ==> dnase-rep-corr
              |                       |                   |                  |                    |
              V                       V                   V                  V                    V
OUTPUTS: tech-pe.bam(c)       filtered-pe.bam(d)  sampled-pe.bam     hs_narrowPeak.bb,bed(e)  corr.txt
                              filtered-pe-qc.txt  sampled-pe-qc.txt  hs_broadPeak.bb,bed.gz
                                                                     hotspot_signal.bw(f)
                                                                     hotspot_qc.txt
```

---------
*Single-end pipeline:*
```
INPUTS:  reads.fq.gz            n*tech-se.bam(g)  filtered-se.bam(h)  filtered-se.bam(h)   2*hotspot_signal.bw(f)
         bwa_index.tgz(a)       n*tech-pe.bam(c)          |           mappability.tgz(b)          |
              |                       |                   |           chrom.sizes                 |
              |                       |                   |                  |                    |
              V                       V                   V                  V                    V
STEPS:   dnase-align-bwa-se ==> bam-filter-se ==> dnase-eval-bam ==> dnase-call-hotspots ==> dnase-rep-corr
              |                       |                   |                  |                    |
              V                       V                   V                  V                    V
OUTPUTS: tech-se.bam(g)       filtered-se.bam(h)  sampled-se.bam     hs_narrowPeak.bb,bed(e)  corr.txt
                              filtered-se-qc.txt  sampled-se-qc.txt  hs_broadPeak.bb,bed.gz
                                                                     hotspot_signal.bw(f)
                                                                     hotspot_qc.txt
```

