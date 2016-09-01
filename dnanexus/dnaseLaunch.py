#!/usr/bin/env python
# dnaseLaunch.py 0.0.1

import argparse,os, sys, json

import dxpy
from launch import Launch
#from template import Launch # (does not use dxencode at all)

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
                "ORDER": { "se": [ "dnase-align-bwa-se" ],
                           "pe": [ "dnase-align-bwa-pe" ] },
                "STEPS": {
                            "dnase-align-bwa-se": {
                                "inputs": { "reads": "reads", "bwa_index": "bwa_index" },
                                "app": "dnase-align-bwa-se", 
                                "params": {"nthreads": "nthreads", "barcode": "barcode" }, 
                                "results": {
                                    "bam_techrep":      "bam_bwa", 
                                    "bam_techrep_qc":   "bam_bwa_qc",
                                },
                                #"output_values": { "all_reads": "all_reads" },
                            },
                            "dnase-align-bwa-pe": {
                                "inputs": { "reads1": "reads1", "reads2": "reads2", "bwa_index": "bwa_index" }, 
                                "app": "dnase-align-bwa-pe", 
                                "params": { "nthreads": "nthreads", "barcode": "barcode", "umi": "umi" }, 
                                "results": {
                                    "bam_techrep":      "bam_bwa", 
                                    "bam_techrep_qc":   "bam_bwa_qc",
                                },
                                #"output_values": { "all_reads": "all_reads" },
                            }, 
                }
        },
        "BIO_REP":  {
                "ORDER": { "se": [  "dnase-filter-se", 
                                    "dnase-eval-bam-alt", 
                                    "dnase-call-hotspots-alt", 
                                 ],
                           "pe": [  "dnase-filter-pe", 
                                    "dnase-eval-bam", 
                                    "dnase-call-hotspots", 
                                 ] },
                "STEPS": {
                            "dnase-filter-pe": {
                                "inputs": { "bam_ABC":    "bam_set" },
                                "app": "dnase-filter-pe", 
                                "params": { "map_thresh": "map_thresh", "nthreads": "nthreads", "umi": "umi" }, 
                                "results": {
                                    "bam_filtered":         "bam_filtered", 
                                    "bam_filtered_qc":      "bam_filtered_qc", 
                                },
                                #"output_values": { "reads_filtered": "filtered_mapped_reads" },
                            },
                            "dnase-filter-se": {
                                "inputs": { "bam_ABC":    "bam_set" },
                                "app": "dnase-filter-se", 
                                "params": { "map_thresh": "map_thresh", "nthreads": "nthreads", "umi": "umi" }, 
                                "results": {
                                    "bam_filtered":         "bam_filtered", 
                                    "bam_filtered_qc":      "bam_filtered_qc", 
                                },
                                #"output_values": { "reads_filtered": "filtered_mapped_reads" },
                            }, 
                            "dnase-eval-bam": {
                                "inputs": { "bam_filtered": "bam_filtered" }, 
                                "app": "dnase-eval-bam", 
                                "params": { "pe_or_se": "pe_or_se", "sample_size": "sample_size", "nthreads": "nthreads" }, 
                                "results": {
                                    "bam_sample":           "bam_sample", 
                                    "bam_sample_qc":        "bam_sample_qc", 
                                },
                                #"output_values": { "sampled_": "sampled_reads" },
                            },
                            "dnase-eval-bam-alt": {
                                "inputs": { "bam_filtered": "bam_filtered" }, 
                                "app": "dnase-eval-bam-se", 
                                "params": { "sample_size": "sample_size" }, 
                                "results": {
                                    #"bam_no_chrM":          "bam_no_chrM", 
                                    #"bam_no_chrM_qc":       "bam_no_chrM_qc", 
                                    "bam_sample":           "bam_sample", 
                                    "bam_sample_qc":        "bam_sample_qc", 
                                },
                                #"output_values": { "sampled_": "sampled_reads" },
                            }, 
                            "dnase-call-hotspots": {
                                "inputs": { "bam_filtered": "bam_to_call", 
                                            "chrom_sizes": "chrom_sizes", 
                                            "hotspot_mappable": "hotspot_mappable" }, 
                                "app": "dnase-call-hotspots", 
                                #"params": {  }, 
                                "results": {
                                    "bed_hotspots":  "bed_hotspots", 
                                     "bb_hotspots":   "bb_hotspots", 
                                    "bed_peaks":     "bed_peaks", 
                                     "bb_peaks":      "bb_peaks",
                                     "bw_density":     "bw_density", 
                                    "hotspots_qc":    "hotspots_qc"
                                },
                            }, 
                            "dnase-call-hotspots-alt": {
                                "inputs": { "bam_filtered": "bam_to_call", "chrom_sizes": "chrom_sizes" }, 
                                "app": "dnase-call-hotspots-alt", 
                                #"params": {  }, 
                                "results": {
                                    "bed_hotspots":  "bed_hotspots", 
                                     "bb_hotspots":   "bb_hotspots", 
                                    "bed_peaks":     "bed_peaks", 
                                     "bb_peaks":      "bb_peaks",
                                     "bw_density":    "bw_density", 
                                    "hotspots_qc":   "hotspots_qc"
                                },
                            }, 
                }
        },
        "COMBINED_REPS": {
                "ORDER": { "se": [ "dnase-rep-corr-alt" ],
                          "pe": [ "dnase-rep-corr",    ] },
                "STEPS": {
                            "dnase-rep-corr": {
                                "inputs": { "density_a":"density_a", "density_b":"density_b" }, 
                                #"app": "dnase-rep-corr", 
                                "results": { "corr_txt": "corr_txt" },
                            },
                            "dnase-rep-corr-alt": {
                                "inputs": { "density_a":"density_a", "density_b":"density_b" }, 
                                #"app": "dnase-rep-corr-alt", 
                                "results": { "corr_txt": "corr_txt" },
                            },
        #                    "dnase-idr": {
        #                        "inputs": {
        #                            "peaks_a":    "peaks_a",    "peaks_b":    "peaks_b", 
        #                            "chrom_sizes": "chrom_sizes" 
        #                        }, 
        #                        "app": "dnase-idr", 
        #                        "params": {}, 
        #                        "results": { "idr_png": "idr_png",
        #                                     "idr_bb":  "idr_bb",
        #                                     "idr_bed": "idr_bed",
        #                        },
        #                    },
        #                    "dnase-idr-alt": {
        #                        "inputs": {
        #                            "peaks_a":    "peaks_a",    "peaks_b":    "peaks_b", 
        #                            "chrom_sizes": "chrom_sizes" 
        #                        }, 
        #                        "app": "dnase-idr-alt", 
        #                        "params": {}, 
        #                        "results": { "idr_png": "idr_png",
        #                                     "idr_bb":  "idr_bb",
        #                                     "idr_bed": "idr_bed",
        #                        },
        #                    },
        #                    "dnase-pool-bioreps": {
        #                        "inputs": {
        #                               "bam_A":    "bam_A",    "bam_B":    "bam_B", 
        #                            "signal_A": "signal_A", "signal_B": "signal_B", 
        #                             "peaks_A":  "peaks_A",  "peaks_B":  "peaks_B", 
        #                            "chrom_sizes": "chrom_sizes" 
        #                        }, 
        #                        "app": "dnase-pool-bioreps", 
        #                        "params": {}, 
        #                        "results": {
        #                            "bam_pooled":       "bam_pooled", 
        #                            "bed_merged":       "bed_merged", 
        #                            "bb_merged":        "bb_merged", 
        #                            "pooled_qc":        "pooled_qc", 
        #                        },
        #                        #"output_values": { "reads_pooled": "reads" },
        #                    },
        #                    "dnase-pool-bioreps-alt": {
        #                        "inputs": {
        #                               "bam_A":    "bam_A",    "bam_B":    "bam_B", 
        #                            "signal_A": "signal_A", "signal_B": "signal_B", 
        #                             "peaks_A":  "peaks_A",  "peaks_B":  "peaks_B", 
        #                            "chrom_sizes": "chrom_sizes" 
        #                        }, 
        #                        "app": "dnase-pool-bioreps-alt", 
        #                        "params": {}, 
        #                        "results": {
        #                            "bam_pooled":       "bam_pooled", 
        #                            "bed_merged":       "bed_merged", 
        #                            "bb_merged":        "bb_merged", 
        #                            "pooled_qc":        "pooled_qc", 
        #                        },
        #                        #"output_values": { "reads_pooled": "reads" },
        #                    },
                }
        }
    }

    FILE_GLOBS = {
        #"reads":                    "/*.fq.gz",
        #"reads1":                   "/*.fq.gz",
        #"reads2":                   "/*.fq.gz",
        # dnase-align-pe/se results:
        "bam_techrep":              "/*_bwa_techrep.bam", 
        "bam_techrep_qc":           "/*_bwa_techrep_qc.txt",
        # dnase-filter-pe/se input/results:
        "bam_ABC":                  "/*_bwa_techrep.bam", 
        "bam_filtered":             "/*_filtered.bam", 
        "bam_filtered_qc":          "/*_filtered_qc.txt", 
        # dnase-eval-bam-pe/se input/results:
        #"bam_filtered": {"bam_filtered": "/*_filtered.bam"}, 
        "bam_sample":               "/*_sample.bam", 
        "bam_sample_qc":            "/*_sample_qc.txt",
        # biorep-call-hotspots input/results:
        #"bam_to_call": {"bam_filtered": "/*_filtered.bam"}, 
        "bed_hotspots":             "/*_hotspots.bed.gz", 
        "bb_hotspots":              "/*_hotspots.bb", 
        "bed_peaks":                "/*_peaks.bed.gz", 
        "bb_peaks":                 "/*_peaks.bb",
        "bw_density":               "/*_density.bw", 
        "hotspots_qc":              "/*_hotspots_qc.txt", 
        # dnase-rep-corr input/results:
        "density_a":                "/*_density.bw",
        "density_b":                "/*_density.bw",
        "corr_txt":                 "/*_corr.txt",
        # dnase-idr input/results:
        #"peaks_a":                  "/*_peaks.bed.gz", # "/*_bwa_biorep_filtered_peaks.bed.gz",
        #"peaks_b":                  "/*_peaks.bed.gz", # "/*_bwa_biorep_filtered_peaks.bed.gz",
        #"idr_bed":                  "/*_idr.bed.gz",
        #"idr_bb":                   "/*_idr.bb",
        #"idr_png":                  "/*_idr.png",
        # dnase-pool-bioreps input/results:
        #"bam_A":                    "/*_filtered.bam",
        #"bam_B":                    "/*_filtered.bam",
        #"signal_A":                 "/*_no_chrM_signal_hotspot.bw",
        #"signal_B":                 "/*_no_chrM_signal_hotspot.bw",
        #"peaks_A":                  "/*_no_chrM_narrowPeak_hotspot.bb",
        #"peaks_B":                  "/*_no_chrM_narrowPeak_hotspot.bb",
        #"bam_pooled":               "/*_pooled.bam", 
        #"bed_merged":               "/*_merged_narrowPeak.bed", 
        #"bb_merged":                "/*_merged_narrowPeak.bb", 
        #"pooled_qc":                "/*_pooled_qc.txt", 
        # pooled-call-hotspots input/results:
        #"bam_to_call":              "/*_pooled.bam", 
        #"pr_bb_hotspot_broadPeak":  "/*_pooled_broadPeak_hotspot.bb", 
        #"pr_bed_hotspot_broadPeak": "/*_pooled_broadPeak_hotspot.bed", 
        #"pr_bed_hotspot_narrowPeak":"/*_pooled_narrowPeak_hotspot.bed", 
        #"pr_bb_hotspot_narrowPeak": "/*_pooled_narrowPeak_hotspot.bb",
        #"pr_bw_hotspot_signal":     "/*_pooled_signal_hotspot.bw", 
        #"pr_bam_hotspot_qc":        "/*_pooled_hotspot_qc.txt",         
    }

    REFERENCE_FILES = {
        # For looking up reference file names.
        # TODO: should use ACCESSION based fileNames
        "bwa_index":   {
                        "GRCh38": "GRCh38_bwa_index.tgz",
                        "hg19":   "hg19_bwa_index.tgz",
                        "mm10":   "mm10_bwa_index.tgz",
                        },
        "hotspot_mappable":   {
                        "GRCh38": "GRCh38_hotspot_mappable.tgz",
                        "hg19":   "hg19_hotspot_mappable.tgz",
                        "mm10":   "mm10_hotspot_mappable.tgz",
                        },
        "chrom_sizes":   {
                        "GRCh38": "GRCh38_EBV.chrom.sizes",
                        "hg19":   "male.hg19.chrom.sizes",
                        "mm10":   "male.mm10.chrom.sizes",
                        }
        }


    def __init__(self):
        Launch.__init__(self)
        self.detect_umi = True     # Only in DNase is there a umi setting buried in the fastq metadata.
        
    def get_args(self):
        '''Parse the input arguments.'''
        ap = Launch.get_args(self,parse=False)
        
        ap.add_argument('-rl', '--read_length',
                        help='The length of reads.',
                        type=int,
                        choices=['32', '36', '40', '50', '58', '72', '76', '100'],
                        default='100',
                        required=False)

        ap.add_argument('--umi',
                        help='Treat fastqs as having UMI flags embedded.',
                        action='store_true',
                        required=False)

        # NOTE: Could override get_args() to have this non-generic control message
        #ap.add_argument('-c', '--control',
        #                help='The control bam for peak calling.',
        #                required=False)

        return ap.parse_args()

    def pipeline_specific_vars(self,args,verbose=False):
        '''Adds pipeline specific variables to a dict, for use building the workflow.'''
        psv = Launch.pipeline_specific_vars(self,args)
        
        # Some specific settings
        psv['nthreads']    = 8
        psv['map_thresh']  = 3
        psv['sample_size'] = 15000000
        psv['read_length'] = args.read_length
        psv['pe_or_se'] = "se"
        if psv['paired_end']:
            psv['pe_or_se'] = "pe"
        if args.umi:
            psv['umi'] = "yes"
        psv['upper_limit'] = 0
        
        if verbose:
            print "Pipeline Specific Vars:"
            print json.dumps(psv,indent=4)
        return psv


    def find_ref_files(self,priors):
        '''Locates all reference files based upon organism and gender.'''
        # TODO:  move all ref files to ref project and replace "/ref/" and self.REF_PROJECT_DEFAULT
        bwa_path = self.psv['refLoc']+"dnase/"+self.REFERENCE_FILES['bwa_index'][self.psv['genome']]
        #base_dir = "/ref/" + self.psv['genome'] + "/dnase/"
        #bwa_path = base_dir+self.REFERENCE_FILES['bwa_index'][self.psv['genome']][self.psv['gender']]
        bwa_fid = self.find_file(bwa_path,self.REF_PROJECT_DEFAULT)
        if bwa_fid == None:
            sys.exit("ERROR: Unable to locate BWA index file '" + self.REF_PROJECT_DEFAULT + ':' + bwa_path + "'")
        else:
            priors['bwa_index'] = bwa_fid

        hot_map = self.psv['refLoc']+"dnase/"+self.REFERENCE_FILES['hotspot_mappable'][self.psv['genome']]
        hot_map_fid = self.find_file(hot_map,self.REF_PROJECT_DEFAULT)
        if hot_map_fid == None:
            sys.exit("ERROR: Unable to locate hotspot_mappable file '" + hot_map + "'")
        else:
            priors['hotspot_mappable'] = hot_map_fid

        chrom_sizes = self.psv['refLoc']+self.REFERENCE_FILES['chrom_sizes'][self.psv['genome']]
        chrom_sizes_fid = self.find_file(chrom_sizes,self.REF_PROJECT_DEFAULT)
        if chrom_sizes_fid == None:
            sys.exit("ERROR: Unable to locate Chrom Sizes file '" + chrom_sizes + "'")
        else:
            priors['chrom_sizes'] = chrom_sizes_fid

        self.psv['ref_files'] = self.REFERENCE_FILES.keys()
        return priors
    

    def add_combining_reps(self, psv):
        '''Defines how replicated are combined.'''
        # OVERRIDING parent because DNase-seq pipeline doesn't follow the standard replicate combination model
        debug=False
        
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
                if debug:
                    print "DEBUG: rep: " + rep_id
                    
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
                        if river.get('umi') == None:  # for DNase: bioreps should combine tech_reps of the same umi setting!
                            river['umi'] = tributary.get('umi')
                        elif river['umi'] != tributary.get('umi'):
                            print >> sys.stderr, "ERROR: mixed UMI setting on tech_reps for bio_rep %s." % river['br']
                            sys.exit(1)
                        if len(river['tributaries']) > 0:
                            river['rep_tech'] += '.'
                        river['rep_tech'] += tributary['rep_tech'][5:]
                        river['tributaries'].append(tributary_id)
            assert len(river['tributaries']) >= 1  # It could be the case that there is one tech_rep for a bio_rep!
            # Try to contract reps1_1 (as opposed to reps1_1.2) to rep1_1
            if '.' not in river['rep_tech']:
                river['rep_tech'] = "rep" + river['rep_tech'][4:]
            # river_id for ['a','b'] = 'b-bio_rep1'
            river_id = river['tributaries'][-1] + '-bio_rep' + str(bio_rep)
            reps[river_id] = river
            # Special case of 2 allows for designating sisters
            if len(river['tributaries']) == 2:
                reps[river['tributaries'][0]]['sister'] = river['tributaries'][1]
                reps[river['tributaries'][1]]['sister'] = river['tributaries'][0]
            if debug:
                print "DEBUG: biorep: " + river_id + " tributaries: " + str(len(river['tributaries']))

        # Finally a pair of bio_reps are merged and processing finishes up
        if 'COMBINED_REPS' in self.PIPELINE_BRANCH_ORDER:
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
                    if '.' not in tributary['rep_tech']:
                        sea['rep_tech'] += tributary['rep_tech'][3:]
                    else:
                        sea['rep_tech'] += tributary['rep_tech'][4:]
                    sea['tributaries'].append(tributary_id)
            
                psv['rep_tech'] = sea['rep_tech']
                reps[self.SEA_ID] = sea
                # Special case of 2 allows for designating sisters
                reps[sea['tributaries'][0]]['sister'] = sea['tributaries'][1]
                reps[sea['tributaries'][1]]['sister'] = sea['tributaries'][0]
                if debug:
                    print "DEBUG: experiment: " + self.SEA_ID + " tributaries: " + str(len(sea['tributaries']))
            elif debug:
                print "DEBUG: Found " + str(len(bio_reps)) + " bio_reps.  If exactly two, they would be combined."
            #print json.dumps(reps,indent=4,sort_keys=True)

    #######################


if __name__ == '__main__':
    '''Run from the command line.'''
    dnaseLaunch = DnaseLaunch()
    dnaseLaunch.run()

