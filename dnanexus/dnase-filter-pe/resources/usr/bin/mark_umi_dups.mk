# From: https://github.com/StamLab/stampipes/makefiles/umi/mark_duplicates.mk
# Used for marking UMI duplicates

###############
# These variables should be passed in for the makefile to work properly
###############
# INPUT_BAM_FILE=
# OUTPUT_BAM_FILE=
# THREADS=
# MAX_MEM=
###############

TMPDIR ?= $(shell pwd)

SAMTOOLS_FILTER_OPTIONS ?= -f 2 -F 512

# Check if we can sort in parallel

SORT := sort

SORT_PARALLEL = $(shell echo `sort --version | grep ^sort | sed 's/^.* //g'` \>= 8.5 | bc)
ifeq ($(SORT_PARALLEL), 1)
	SORT += --parallel=$(THREADS)
endif

# Allows us to split up the workload
CHROMS ?= $(shell samtools view -H $(INPUT_BAM_FILE) | grep '^@SQ'| sed 's/^@SQ\s\+SN:\([^\s]\+\)\s\+.*/\1/' )
CHROM_FILES = $(addprefix $(TMPDIR)/,$(addsuffix .marked.bam,$(CHROMS)))

INDEX ?= $(INPUT_BAM_FILE).bai

all: markdups

markdups: $(OUTPUT_BAM_FILE)

$(OUTPUT_BAM_FILE): $(CHROM_FILES) $(TMPDIR)/unmapped.marked.bam
	samtools merge -@ $(THREADS) $@ $^

$(TMPDIR)/%.marked.bam: $(TMPDIR)/%.reads.bam $(TMPDIR)/%.reads.dups
	( \
		samtools view -H $(TMPDIR)/$*.reads.bam; \
	\
		samtools view $(TMPDIR)/$*.reads.bam \
			| join -j 1 -v 1 - $(TMPDIR)/$*.reads.dups \
			| tr " " "\t" \
			| awk -v OFS="\t" '{ $$2 = and($$2, compl(lshift(1, 10))); print; }' ;\
	\
		samtools view $(TMPDIR)/$*.reads.bam \
			| join -j 1 - $(TMPDIR)/$*.reads.dups \
			| tr " " "\t" \
			| awk -v OFS="\t" '{ $$2 = or($$2, lshift(1, 10)); print; }' ;\
	) \
		| samtools sort -m $(MAX_MEM)G -@ $(THREADS) -O bam -T $(TMPDIR)/$*.tmpsort -o $@ -
	rm $^ $(TMPDIR)/$*.reads.sorted $(TMPDIR)/$*.reads.alldups $(TMPDIR)/$*.reads.firstdup

# Prints out all duplicates (except the first one) (previous step
# sorts reads by highest MAPQ
#
# Strategy is to find all the duplicate lines, then remove the first
# instance of a duplicate, and return the rest

$(TMPDIR)/%.reads.dups: $(TMPDIR)/%.reads.alldups $(TMPDIR)/%.reads.firstdup
	join -j 1 -v 1 $(TMPDIR)/$*.reads.alldups $(TMPDIR)/$*.reads.firstdup > $@

$(TMPDIR)/%.reads.firstdup: $(TMPDIR)/%.reads.sorted
	cat $^ | uniq -f2 -d | cut -f1 | $(SORT) --buffer-size=$(MAX_MEM)G > $@

$(TMPDIR)/%.reads.alldups: $(TMPDIR)/%.reads.sorted
	cat $^ | uniq -f2 -D | cut -f1 | $(SORT) --buffer-size=$(MAX_MEM)G > $@

# Sort the fragments by position, strand, UMI and then quality
# TAG ID, MAPQ, chr, start, end, strand, UMI

$(TMPDIR)/%.reads.sorted: $(TMPDIR)/%.reads.bam
	samtools view $(SAMTOOLS_FILTER_OPTIONS) $^ \
		| awk -f /usr/bin/umi_sort_sam_annotate.awk \
		| $(SORT) --buffer-size=$(MAX_MEM)G -k3,3 -k4,4g -k5,5g -k6,6 -k7,7 -k2,2gr \
	>$@

# Sort by name to start
# Splits by sequence chr1,chr2,etc
$(TMPDIR)/%.reads.bam: $(INPUT_BAM_FILE) $(INDEX)
	( samtools view -H $<  ;\
		samtools view $< $* | $(SORT) --buffer-size=$(MAX_MEM)G -k1,1 ; \
	) \
	| samtools view -S -1 - -o $@

$(TMPDIR)/unmapped.marked.bam: $(INPUT_BAM_FILE) $(INDEX)
	samtools view -o $@ -h $< '*'

$(INDEX): $(INPUT_BAM_FILE)
	samtools index $<

