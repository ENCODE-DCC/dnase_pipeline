#!/usr/bin/env python
# dnaseLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
#from dxencode import dxencode as dxencode
import dxencode as dxencode
from launch import Launch

class DnaseLaunch(Launch):
    '''Descendent from Launch class with 'rampage' methods'''

    PIPELINE_NAME = "dnase-seq"
    ''' This must match the assay type returned by dxencode.get_assay_type({exp_id}).'''
    PIPELINE_HELP = "Launches '"+PIPELINE_NAME+"' pipeline " + \
                    "analysis for one replicate or combined replicates. "
    ''' This help title should name pipline and whether combined replicates are supported.'''
                    
    RESULT_FOLDER_DEFAULT = '/dnase/'
    ''' This the default location to place results folders for each experiment.'''
    
    PIPELINE_BRANCH_ORDER = [ "TECH_REP", "BIO_REP", "COMBINED_REPS" ]
    '''A pipeline is frequently made of branches that flow into each other, such as replicate level to combined replicate.'''
    
    PIPELINE_BRANCHES = {
    #'''Each branch must define the 'steps' and their (artificially) linear order.'''
        "TECH_REP": {
                "ORDER": { "se": [ "align-bwa-se" ],
                           "pe": [ "align-bwa-pe" ] },
                "STEPS": {
                            "align-bwa-se": {
                                "inputs": { "reads": "reads", "bwa_index": "bwa_index" },
                                "app": "align-bwa-se", 
                                "params": {"nthreads": "nthreads"}, 
                                "results": {
                                    "bam_bwa":      "bam_bwa", 
                                    "bam_bwa_qc":   "bam_bwa_qc"
                                }
                            },
                            "align-bwa-pe": {
                                "inputs": { "reads1": "reads1", "reads2": "reads2", "bwa_index": "bwa_index" }, 
                                "app": "align-bwa-pe", 
                                "params": { "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_bwa":      "bam_bwa", 
                                    "bam_bwa_qc":   "bam_bwa_qc"
                                }
                            }, 
                }
        },
        "BIO_REP":  {
                "ORDER": { "se": [ "merge-bams", "bam-filter-se", "call-hotspots", "sample-hotspots" ],
                           "pe": [ "merge-bams", "bam-filter-pe", "call-hotspots", "sample-hotspots" ] },
                "STEPS": {
                            "merge-bams": {
                                "inputs": { "bam_set":    "bam_set" },
                                "app": "merge-bams", 
                                "params": { "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_merged":      "bam_merged", 
                                    "bam_merged_qc":   "bam_merged_qc"
                                }
                            },
                            "bam-filter-pe": {
                                "inputs": { "bam_merged": "bam_bwa" }, 
                                "app": "bam-filter-pe", 
                                "params": { "sample_size": "sample_size", "map_thresh": "map_thresh", "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_filtered_qc":      "bam_filtered_qc", 
                                    "bam_sample":           "bam_sample", 
                                    "bam_filtered_qc_full": "bam_filtered_qc_full", 
                                    "bam_filtered":         "bam_filtered", 
                                    "bam_sample_pbc":       "bam_sample_pbc", 
                                    "bam_no_chrM":          "bam_no_chrM", 
                                    "bam_sample_stats":     "bam_sample_stats", 
                                    "bam_sample_spp":       "bam_sample_spp"
                                }
                            },
                            "bam-filter-se": {
                                "inputs": { "bam_merged": "bam_bwa" }, 
                                "app": "bam-filter-se", 
                                "params": { "sample_size": "sample_size", "map_thresh": "map_thresh", "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_filtered_qc":      "bam_filtered_qc", 
                                    "bam_sample":           "bam_sample", 
                                    "bam_filtered_qc_full": "bam_filtered_qc_full", 
                                    "bam_filtered":         "bam_filtered", 
                                    "bam_sample_pbc":       "bam_sample_pbc", 
                                    "bam_no_chrM":          "bam_no_chrM", 
                                    "bam_sample_stats":     "bam_sample_stats", 
                                    "bam_sample_spp":       "bam_sample_spp"
                                }
                            }, 
                            "call-hotspots": {
                                "inputs": { "bam_no_chrM": "bam_to_call", "chrom_sizes": "chrom_sizes" }, 
                                "app": "call-hotspots", 
                                "params": { "read_length": "read_length", "genome": "genome" }, 
                                "results": {
                                     "rs_bb_hotspot_broadPeak":   "bb_hotspot_broadPeak", 
                                    "rs_bed_hotspot_broadPeak":  "bed_hotspot_broadPeak", 
                                     "rs_bb_hotspot_narrowPeak":  "bb_hotspot_narrowPeak",
                                    "rs_bed_hotspot_narrowPeak": "bed_hotspot_narrowPeak", 
                                     "rs_bw_hotspot_signal":      "bw_hotspot_signal", 
                                    "rs_bam_hotspot_qc":         "bam_hotspot_qc"
                                }
                            }, 
                            "sample-hotspots": {
                                "inputs": { "bam_no_chrM": "bam_to_sample", "chrom_sizes": "chrom_sizes" }, 
                                "app": "sample-hotspots", 
                                "params": { "read_length": "read_length", "genome": "genome"}, 
                                "results": {
                                    "hotspot_sample_5M_qc": "hotspot_sample_5M_qc", 
                                    "bam_sample_5M":        "bam_sample_5M"
                                }
                            }
                }
        },
        "COMBINED_REPS": {
                "ORDER": [ "merge-replicates", "call-merged-hotspots" ],
                "STEPS": {
                            "merge-replicates": {
                                "inputs": {
                                       "bam_A":    "bam_A",    "bam_B":    "bam_B", 
                                    "signal_A": "signal_A", "signal_B": "signal_B", 
                                     "peaks_A":  "peaks_A",  "peaks_B":  "peaks_B", 
                                    "chrom_sizes": "chrom_sizes" 
                                }, 
                                "app": "merge-replicates", 
                                "params": {}, 
                                "results": {
                                    "bed_merged":       "bed_merged", 
                                    "bam_pooled":       "bam_pooled", 
                                    "signal_corr_qc":   "signal_corr_qc", 
                                    "bb_merged":        "bb_merged", 
                                    "peaks_overlap_qc": "peaks_overlap_qc"
                                }
                            },
                            "call-merged-hotspots": {
                                "inputs": { "bam_pooled": "bam_to_call", "chrom_sizes": "chrom_sizes" }, 
                                "app": "call-hotspots", 
                                "params": { "read_length": "read_length", "genome": "genome" }, 
                                "results": {
                                     "cs_bb_hotspot_broadPeak":   "bb_hotspot_broadPeak", 
                                    "cs_bed_hotspot_broadPeak":  "bed_hotspot_broadPeak", 
                                     "cs_bb_hotspot_narrowPeak":  "bb_hotspot_narrowPeak",
                                    "cs_bed_hotspot_narrowPeak": "bed_hotspot_narrowPeak", 
                                     "cs_bw_hotspot_signal":      "bw_hotspot_signal", 
                                    "cs_bam_hotspot_qc":         "bam_hotspot_qc" 
                                }
                            } 
                }
        }
    }

    # TODO: Remove when we are certain we are abandoning standard replicate combination model
    #PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
    #'''A pipeline is frequently made of branches that flow into each other, such as replicate level to combined replicate.'''
    #
    #PIPELINE_BRANCHES = {
    ##'''Each branch must define the 'steps' and their (artificially) linear order.'''
    #    "REP": {
    #            "ORDER": { "se": [ "align-bwa-se", "bam-filter-se", "call-hotspots", "sample-hotspots" ],
    #                       "pe": [ "align-bwa-pe", "bam-filter-pe", "call-hotspots", "sample-hotspots" ] },
    #            "STEPS": {
    #                        "align-bwa-se": {
    #                            "inputs": { "reads": "reads", "bwa_index": "bwa_index" },
    #                            "app": "align-bwa-se", 
    #                            "params": {"nthreads": "nthreads"}, 
    #                            "results": {
    #                                "bam_bwa":      "bam_bwa", 
    #                                "bam_bwa_qc":   "bam_bwa_qc"
    #                            }
    #                        },
    #                        "align-bwa-pe": {
    #                            "inputs": { "reads1": "reads1", "reads2": "reads2", "bwa_index": "bwa_index" }, 
    #                            "app": "align-bwa-pe", 
    #                            "params": { "nthreads": "nthreads" }, 
    #                            "results": {
    #                                "bam_bwa":      "bam_bwa", 
    #                                "bam_bwa_qc":   "bam_bwa_qc"
    #                            }
    #                        }, 
    #                        "bam-filter-pe": {
    #                            "inputs": { "bam_bwa": "bam_bwa" }, 
    #                            "app": "bam-filter-pe", 
    #                            "params": { "sample_size": "sample_size", "map_thresh": "map_thresh", "nthreads": "nthreads" }, 
    #                            "results": {
    #                                "bam_filtered_qc":      "bam_filtered_qc", 
    #                                "bam_sample":           "bam_sample", 
    #                                "bam_filtered_qc_full": "bam_filtered_qc_full", 
    #                                "bam_filtered":         "bam_filtered", 
    #                                "bam_sample_pbc":       "bam_sample_pbc", 
    #                                "bam_no_chrM":          "bam_no_chrM", 
    #                                "bam_sample_stats":     "bam_sample_stats", 
    #                                "bam_sample_spp":       "bam_sample_spp"
    #                            }
    #                        },
    #                        "bam-filter-se": {
    #                            "inputs": { "bam_bwa": "bam_bwa" }, 
    #                            "app": "bam-filter-se", 
    #                            "params": { "sample_size": "sample_size", "map_thresh": "map_thresh", "nthreads": "nthreads" }, 
    #                            "results": {
    #                                "bam_filtered_qc":      "bam_filtered_qc", 
    #                                "bam_sample":           "bam_sample", 
    #                                "bam_filtered_qc_full": "bam_filtered_qc_full", 
    #                                "bam_filtered":         "bam_filtered", 
    #                                "bam_sample_pbc":       "bam_sample_pbc", 
    #                                "bam_no_chrM":          "bam_no_chrM", 
    #                                "bam_sample_stats":     "bam_sample_stats", 
    #                                "bam_sample_spp":       "bam_sample_spp"
    #                            }
    #                        }, 
    #                        "call-hotspots": {
    #                            "inputs": { "bam_no_chrM": "bam_to_call", "chrom_sizes": "chrom_sizes" }, 
    #                            "app": "call-hotspots", 
    #                            "params": { "read_length": "read_length", "genome": "genome" }, 
    #                            "results": {
    #                                 "rs_bb_hotspot_broadPeak":   "bb_hotspot_broadPeak", 
    #                                "rs_bed_hotspot_broadPeak":  "bed_hotspot_broadPeak", 
    #                                 "rs_bb_hotspot_narrowPeak":  "bb_hotspot_narrowPeak",
    #                                "rs_bed_hotspot_narrowPeak": "bed_hotspot_narrowPeak", 
    #                                 "rs_bw_hotspot_signal":      "bw_hotspot_signal", 
    #                                "rs_bam_hotspot_qc":         "bam_hotspot_qc"
    #                            }
    #                        }, 
    #                        "sample-hotspots": {
    #                            "inputs": { "bam_no_chrM": "bam_to_sample", "chrom_sizes": "chrom_sizes" }, 
    #                            "app": "sample-hotspots", 
    #                            "params": { "read_length": "read_length", "genome": "genome"}, 
    #                            "results": {
    #                                "hotspot_sample_5M_qc": "hotspot_sample_5M_qc", 
    #                                "bam_sample_5M":        "bam_sample_5M"
    #                            }
    #                        }
    #            }
    #    },
    #    "COMBINED_REPS": {
    #            "ORDER": [ "merge-replicates", "call-merged-hotspots" ],
    #            "STEPS": {
    #                        "merge-replicates": {
    #                            "inputs": {
    #                                   "bam_A":    "bam_A",    "bam_B":    "bam_B", 
    #                                "signal_A": "signal_A", "signal_B": "signal_B", 
    #                                 "peaks_A":  "peaks_A",  "peaks_B":  "peaks_B", 
    #                                "chrom_sizes": "chrom_sizes" 
    #                            }, 
    #                            "app": "merge-replicates", 
    #                            "params": {}, 
    #                            "results": {
    #                                "bed_merged":       "bed_merged", 
    #                                "bam_pooled":       "bam_pooled", 
    #                                "signal_corr_qc":   "signal_corr_qc", 
    #                                "bb_merged":        "bb_merged", 
    #                                "peaks_overlap_qc": "peaks_overlap_qc"
    #                            }
    #                        },
    #                        "call-merged-hotspots": {
    #                            "inputs": { "bam_pooled": "bam_to_call", "chrom_sizes": "chrom_sizes" }, 
    #                            "app": "call-hotspots", 
    #                            "params": { "read_length": "read_length", "genome": "genome" }, 
    #                            "results": {
    #                                 "cs_bb_hotspot_broadPeak":   "bb_hotspot_broadPeak", 
    #                                "cs_bed_hotspot_broadPeak":  "bed_hotspot_broadPeak", 
    #                                 "cs_bb_hotspot_narrowPeak":  "bb_hotspot_narrowPeak",
    #                                "cs_bed_hotspot_narrowPeak": "bed_hotspot_narrowPeak", 
    #                                 "cs_bw_hotspot_signal":      "bw_hotspot_signal", 
    #                                "cs_bam_hotspot_qc":         "bam_hotspot_qc" 
    #                            }
    #                        } 
    #            }
    #    }
    #}

    FILE_GLOBS = {
        #"reads":                    "/*.fq.gz",
        #"reads1":                   "/*.fq.gz",
        #"reads2":                   "/*.fq.gz",
        "bam_bwa":                  "/*_bwa.bam", 
        "bam_bwa_qc":               "/*_bam_qc.txt",
        "bam_filtered":             "/*_filtered.bam", 
        "bam_filtered_qc":          "/*_filtered_qc.txt", 
        "bam_filtered_qc_full":     "/*_filtered_qc_full.txt", 
        "bam_merged":               "/*_merged.bam", 
        "bam_merged_qc":            "/*_merged_bam_qc.txt", 
        "bam_no_chrM":              "/*_no_chrM.bam", 
        "bam_sample":               "/*_sample.bam", 
        "bam_sample_pbc":           "/*_sample_pbc.txt", 
        "bam_sample_spp":           "/*_sample_spp.txt", 
        "bam_sample_stats":         "/*_sample_stats.txt", 
        "bam_sample_5M":            "/*_sample_5M.bam", 
        "bam_set":                  "/*_bwa.bam", 
        "bam_pooled":               "/*_pooled.bam", 
        #"bam_to_call":              "/*_pooled.bam", 
        "rs_bb_hotspot_broadPeak":  "/*_no_chrM_broadPeak_hotspot.bb", 
        "rs_bed_hotspot_broadPeak": "/*_no_chrM_broadPeak_hotspot.bed", 
        "rs_bed_hotspot_narrowPeak":"/*_no_chrM_narrowPeak_hotspot.bed", 
        "rs_bb_hotspot_narrowPeak": "/*_no_chrM_narrowPeak_hotspot.bb",
        "rs_bw_hotspot_signal":     "/*_no_chrM_signal_hotspot.bw", 
        "rs_bam_hotspot_qc":        "/*_no_chrM_hotspot_qc.txt", 
        "cs_bb_hotspot_broadPeak":  "/*_pooled_broadPeak_hotspot.bb", 
        "cs_bed_hotspot_broadPeak": "/*_pooled_broadPeak_hotspot.bed", 
        "cs_bed_hotspot_narrowPeak":"/*_pooled_narrowPeak_hotspot.bed", 
        "cs_bb_hotspot_narrowPeak": "/*_pooled_narrowPeak_hotspot.bb",
        "cs_bw_hotspot_signal":     "/*_pooled_signal_hotspot.bw", 
        "cs_bam_hotspot_qc":        "/*_pooled_hotspot_qc.txt", 
        "bb_merged":                "/*_narrowPeak_merged.bb", 
        "bed_merged":               "/*_narrowPeak_merged.bed", 
        "hotspot_sample_5M_qc":     "/*_sample_5M_hotspot_qc.txt", 
        "signal_corr_qc":           "/*_signal_corr_qc.txt", 
        "peaks_overlap_qc":         "/*_narrowPeaks_overlap_qc.txt",
        "bam_A":                    "/*_filtered.bam",
        "bam_B":                    "/*_filtered.bam",
        "signal_A":                 "/*_no_chrM_signal_hotspot.bw",
        "signal_B":                 "/*_no_chrM_signal_hotspot.bw",
        "peaks_A":                  "/*_no_chrM_narrowPeak_hotspot.bb",
        "peaks_B":                  "/*_no_chrM_narrowPeak_hotspot.bb"
    }

    REF_PROJECT_DEFAULT = "scratchPad"  # TODO: move all ref files to ref project!
    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should use ACCESSION based fileNames
        "bwa_index":   {
                        "hg19": {
                                "female":   "hg19_female_bwa_index.tgz",
                                "male":     "hg19_male_bwa_index.tgz"
                                }
                        },
        "chrom_sizes":   {
                        "hg19": {
                                "female":   "female.hg19.chrom.sizes",
                                "male":     "male.hg19.chrom.sizes"
                                },
                        "mm10": {
                                "female":   "female.mm10.chrom.sizes",
                                "male":     "male.mm10.chrom.sizes"
                                }
                        }
        }


    def __init__(self):
        Launch.__init__(self)
        
    def get_args(self):
        '''Parse the input arguments.'''
        ap = Launch.get_args(self,parse=False)
        
        ap.add_argument('-rl', '--read_length',
                        help='The length of reads.',
                        type=int,
                        choices=['32', '36', '40', '50', '58', '72', '76', '100'],
                        default='100',
                        required=False)

        # NOTE: Could override get_args() to have this non-generic control message
        #ap.add_argument('-c', '--control',
        #                help='The control bam for peak calling.',
        #                required=False)

        return ap.parse_args()

    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        psv = Launch.pipeline_specific_vars(self,args)
        
        #if not psv['paired_end']:
        #    print "Rampage is always expected to be paired-end but mapping says otherwise."
        #    sys.exit(1)

        # Some specific settings
        psv['nthreads']   = 8
        psv['map_thresh']   = 3
        psv['sample_size']   = 15000000
        psv['read_length']   = args.read_length
        
        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        # TODO:  move all ref files to ref project and replace "/ref/" and self.REF_PROJECT_DEFAULT
        #bwaIx = self.psv['refLoc']+self.REFERENCE_FILES['bwa_index'][self.psv['genome']][self.psv['gender']]
        bwaIx = "/ref/"+self.REFERENCE_FILES['bwa_index'][self.psv['genome']][self.psv['gender']]
        #bwaIxFid = dxencode.find_file(bwaIx,dxencode.REF_PROJECT_DEFAULT)
        bwaIxFid = dxencode.find_file(bwaIx,self.REF_PROJECT_DEFAULT)
        if bwaIxFid == None:
            sys.exit("ERROR: Unable to locate BWA index file '" + bwaIx + "'")
        else:
            priors['bwa_index'] = bwaIxFid

        chromSizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']][self.psv['gender']]
        chromSizesFid = dxencode.find_file(chromSizes,dxencode.REF_PROJECT_DEFAULT)
        if chromSizesFid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chromSizes + "'")
        else:
            priors['chrom_sizes'] = chromSizesFid
        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
    

    def add_combining_reps(self, psv):
        '''Defines how replicated are combined.'''
        # OVERRIDING parent because DNase-seq pipeline doesn't follow the standard replicate combination model
        
        reps = psv['reps']
        # In the 'standard combining model' PIPELINE_BRANCH_ORDER = [ "REP", "COMBINED_REPS" ]
        # and all replicates are in psv['reps'] keyed as 'a','b',... and having rep['rep_tech'] = 'rep1_1'
        # All these simple reps will have rep['branch_id'] = "REP"
        
        # First, each tech_rep is processed individually
        bio_reps = []
        for rep_id in sorted( reps.keys() ):
            if len(rep_id) == 1: # single letter: simple replicate
                rep = reps[rep_id]
                rep['branch_id'] = "TECH_REP"
                if rep['br'] not in bio_reps:
                    bio_reps.append(rep['br'])
                else:
                    self.combined_reps = True  # More than one tech_rep per bio_rep so combining will be done!
                    
        # Next bio_reps have their technical replicates merged and processing continues
        for bio_rep in bio_reps:
            river = {}
            river['branch_id'] = "BIO_REP"
            river['tributaries'] = []
            river['rep_tech'] = 'reps' + str(bio_rep) + '_'  # reps1_1.2.3 is rep1_1 + rep1_2 + rep1_3
            river['br'] = bio_rep
            for tributary_id in sorted( reps.keys() ): 
                if len(tributary_id) == 1:
                    tributary = reps[tributary_id]
                    if tributary['br'] == bio_rep:
                        if len(river['tributaries']) > 0:
                            river['rep_tech'] += '.'
                        river['rep_tech'] += tributary['rep_tech'][5:]
                        river['tributaries'].append(tributary_id)
            assert len(river['tributaries']) >= 1  # It could be the case that there is one tech_rep for a bio_rep!
            # river_id for ['a','b'] = 'b-bio_rep1'
            river_id = river['tributaries'][-1] + '-bio_rep' + str(bio_rep)
            reps[river_id] = river

        # Finally a pair of bio_reps are merged and processing finishes up
        if len(bio_reps) == 2:
            self.combined_reps = True  # More than one bio_rep so combining will be done!
            sea = {} # SEA is the final branch into which all tributaries flow
            sea['branch_id'] = 'COMBINED_REPS'
            sea['tributaries'] = []
            sea['rep_tech'] = 'reps'
            for tributary_id in sorted( reps.keys() ):
                if len(tributary_id) == 1:  # ignore the simple reps
                    continue 
                tributary = reps[tributary_id]
                if len(sea['tributaries']) > 0:
                    sea['rep_tech'] += '-'
                sea['rep_tech'] += tributary['rep_tech'][4:]
                sea['tributaries'].append(tributary_id)
        
            psv['rep_tech'] = sea['rep_tech']
            reps[self.SEA_ID] = sea
        #else:
        #    print "Found " + str(len(bio_reps)) + " bio_reps.  If exactly two, they would be combined."
        #print json.dumps(reps,indent=4,sort_keys=True)
            

    #######################


if __name__ == '__main__':
    '''Run from the command line.'''
    dnaseLaunch = DnaseLaunch()
    dnaseLaunch.run()

