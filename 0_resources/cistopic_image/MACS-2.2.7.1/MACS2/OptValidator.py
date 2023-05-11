# Time-stamp: <2019-09-25 10:15:20 taoliu>

"""Module Description

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------
import sys
import os
import re
import logging
from argparse import ArgumentError
from subprocess import Popen, PIPE
from math import log
from MACS2.IO.Parser import BEDParser, ELANDResultParser, ELANDMultiParser, \
    ELANDExportParser, SAMParser, BAMParser, BAMPEParser,\
    BEDPEParser, BowtieParser,  guess_parser
# ------------------------------------
# constants
# ------------------------------------

efgsize = {"hs":2.7e9,
           "mm":1.87e9,
           "ce":9e7,
           "dm":1.2e8}

# ------------------------------------
# Misc functions
# ------------------------------------
def opt_validate ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logging.error("Error when interpreting --gsize option: %s" % options.gsize)
            logging.error("Available shortcuts of effective genome sizes are %s" % ",".join(list(efgsize.keys())))
            sys.exit(1)

    # format
    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
        options.nomodel = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
        options.nomodel = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # duplicate reads
    if options.keepduplicates != "auto" and options.keepduplicates != "all":
        if not options.keepduplicates.isdigit():
            logging.error("--keep-dup should be 'auto', 'all' or an integer!")
            sys.exit(1)

    # shiftsize>0
    #if options.shiftsize:               # only if --shiftsize is set, it's true
    #    options.extsize = 2 * options.shiftsize
    #else:                               # if --shiftsize is not set
    #    options.shiftsize = options.extsize / 2
    if options.extsize < 1 :
        logging.error("--extsize must >= 1!")
        sys.exit(1)

    # refine_peaks, call_summits can't be combined with --broad
    #if options.broad and (options.refine_peaks or options.call_summits):
    #    logging.error("--broad can't be combined with --refine-peaks or --call-summits!")
    #    sys.exit(1)

    if options.broad and options.call_summits:
        logging.error("--broad can't be combined with --call-summits!")
        sys.exit(1)

    if options.pvalue:
        # if set, ignore qvalue cutoff
        options.log_qvalue = None
        options.log_pvalue = log(options.pvalue,10)*-1
    else:
        options.log_qvalue = log(options.qvalue,10)*-1
        options.log_pvalue = None
    if options.broad:
        options.log_broadcutoff = log(options.broadcutoff,10)*-1
    
    # uppercase the format string 
    options.format = options.format.upper()

    # d_min is non-negative
    if options.d_min < 0:
        logging.error("Minimum fragment size shouldn't be negative!" % options.d_min)
        sys.exit(1)
    
    # upper and lower mfold
    options.lmfold = options.mfold[0]
    options.umfold = options.mfold[1]
    if options.lmfold > options.umfold:
        logging.error("Upper limit of mfold should be greater than lower limit!" % options.mfold)
        sys.exit(1)
    
    # output filenames
    options.peakxls = os.path.join( options.outdir, options.name+"_peaks.xls" )
    options.peakbed = os.path.join( options.outdir, options.name+"_peaks.bed" )
    options.peakNarrowPeak = os.path.join( options.outdir, options.name+"_peaks.narrowPeak" )
    options.peakBroadPeak = os.path.join( options.outdir, options.name+"_peaks.broadPeak" )
    options.peakGappedPeak = os.path.join( options.outdir, options.name+"_peaks.gappedPeak" )
    options.summitbed = os.path.join( options.outdir, options.name+"_summits.bed" )
    options.bdg_treat = os.path.join( options.outdir, options.name+"_treat_pileup.bdg" )
    options.bdg_control= os.path.join( options.outdir, options.name+"_control_lambda.bdg" )
    if options.cutoff_analysis:
        options.cutoff_analysis_file = os.path.join( options.outdir, options.name+"_cutoff_analysis.txt" )
    else:
        options.cutoff_analysis_file = "None"
    #options.negxls  = os.path.join( options.name+"_negative_peaks.xls" )
    #options.diagxls = os.path.join( options.name+"_diag.xls" )
    options.modelR  = os.path.join( options.outdir, options.name+"_model.r" )
    #options.pqtable  = os.path.join( options.outdir, options.name+"_pq_table.txt" )

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    options.argtxt = "\n".join((
        "# Command line: %s" % " ".join(sys.argv[1:]),\
        "# ARGUMENTS LIST:",\
        "# name = %s" % (options.name),\
        "# format = %s" % (options.format),\
        "# ChIP-seq file = %s" % (options.tfile),\
        "# control file = %s" % (options.cfile),\
        "# effective genome size = %.2e" % (options.gsize),\
        #"# tag size = %d" % (options.tsize),\
        "# band width = %d" % (options.bw),\
        "# model fold = %s\n" % (options.mfold),\
        ))

    if options.pvalue:
        if options.broad:
            options.argtxt +=  "# pvalue cutoff for narrow/strong regions = %.2e\n" % (options.pvalue)
            options.argtxt +=  "# pvalue cutoff for broad/weak regions = %.2e\n" % (options.broadcutoff)
            options.argtxt +=  "# qvalue will not be calculated and reported as -1 in the final output.\n"
        else:
            options.argtxt +=  "# pvalue cutoff = %.2e\n" % (options.pvalue)
            options.argtxt +=  "# qvalue will not be calculated and reported as -1 in the final output.\n"
    else:
        if options.broad:
            options.argtxt +=  "# qvalue cutoff for narrow/strong regions = %.2e\n" % (options.qvalue)
            options.argtxt +=  "# qvalue cutoff for broad/weak regions = %.2e\n" % (options.broadcutoff)
        else:
            options.argtxt +=  "# qvalue cutoff = %.2e\n" % (options.qvalue)

    if options.maxgap:
        options.argtxt += "# The maximum gap between significant sites = %d\n" % options.maxgap
    else:
        options.argtxt += "# The maximum gap between significant sites is assigned as the read length/tag size.\n" 
    if options.minlen:
        options.argtxt += "# The minimum length of peaks = %d\n" % options.minlen
    else:
        options.argtxt += "# The minimum length of peaks is assigned as the predicted fragment length \"d\".\n"
        
    if options.downsample:
        options.argtxt += "# Larger dataset will be randomly sampled towards smaller dataset.\n"
        if options.seed >= 0:
            options.argtxt += "# Random seed has been set as: %d\n" % options.seed
    else:
        if options.scaleto == "large":
            options.argtxt += "# Smaller dataset will be scaled towards larger dataset.\n"
        else:
            options.argtxt += "# Larger dataset will be scaled towards smaller dataset.\n"

    if options.ratio != 1.0:
        options.argtxt += "# Using a custom scaling factor: %.2e\n" % (options.ratio)
	
    if options.cfile:
        options.argtxt += "# Range for calculating regional lambda is: %d bps and %d bps\n" % (options.smalllocal,options.largelocal)
    else:
        options.argtxt += "# Range for calculating regional lambda is: %d bps\n" % (options.largelocal)

    if options.broad:
        options.argtxt += "# Broad region calling is on\n"
    else:
        options.argtxt += "# Broad region calling is off\n"

    if options.fecutoff != 1.0:
        options.argtxt += "# Additional cutoff on fold-enrichment is: %.2f\n" % (options.fecutoff)

    if options.format in ["BAMPE", "BEDPE"]:
        # neutralize SHIFT
        options.shift = 0
        options.argtxt += "# Paired-End mode is on\n"
    else:
        options.argtxt += "# Paired-End mode is off\n"

    #if options.refine_peaks:
    #    options.argtxt += "# Refining peak for read balance is on\n"
    if options.call_summits:
        options.argtxt += "# Searching for subpeak summits is on\n"

    if options.do_SPMR and options.store_bdg:
        options.argtxt += "# MACS will save fragment pileup signal per million reads\n"        

    return options

def diff_opt_validate ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # format
    options.gzip_flag = False           # if the input is gzip file
    
#    options.format = options.format.upper()
    # fox this stuff
#    if True: pass
#    elif options.format == "AUTO":
#        options.parser = guess_parser
#    else:
#        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
#        sys.exit(1)
    
    if options.peaks_pvalue:
        # if set, ignore qvalue cutoff
        options.peaks_log_qvalue = None
        options.peaks_log_pvalue = log(options.peaks_pvalue,10)*-1
        options.track_score_method = 'p'
    else:
        options.peaks_log_qvalue = log(options.peaks_qvalue,10)*-1
        options.peaks_log_pvalue = None
        options.track_score_method = 'q'

    if options.diff_pvalue:
        # if set, ignore qvalue cutoff
        options.log_qvalue = None
        options.log_pvalue = log(options.diff_pvalue,10)*-1
        options.score_method = 'p'
    else:
        options.log_qvalue = log(options.diff_qvalue,10)*-1
        options.log_pvalue = None
        options.score_method = 'q'
    
    # output filenames
    options.peakxls = options.name+"_diffpeaks.xls"
    options.peakbed = options.name+"_diffpeaks.bed"
    options.peak1xls = options.name+"_diffpeaks_by_peaks1.xls"
    options.peak2xls = options.name+"_diffpeaks_by_peaks2.xls"
    options.bdglogLR = options.name+"_logLR.bdg"
    options.bdgpvalue = options.name+"_logLR.bdg"
    options.bdglogFC = options.name+"_logLR.bdg"
    
    options.call_peaks = True
    if not (options.peaks1 == '' or options.peaks2 == ''):
        if options.peaks1 == '':
            raise ArgumentError('peaks1', 'Must specify both peaks1 and peaks2, or neither (to call peaks again)')
        elif options.peaks2 == '':
            raise ArgumentError('peaks2', 'Must specify both peaks1 and peaks2, or neither (to call peaks again)')
        options.call_peaks = False
        options.argtxt = "\n".join((
            "# ARGUMENTS LIST:",\
            "# name = %s" % (options.name),\
#            "# format = %s" % (options.format),\
            "# ChIP-seq file 1 = %s" % (options.t1bdg),\
            "# control file 1 = %s" % (options.c1bdg),\
            "# ChIP-seq file 2 = %s" % (options.t2bdg),\
            "# control file 2 = %s" % (options.c2bdg),\
            "# Peaks, condition 1 = %s" % (options.peaks1),\
            "# Peaks, condition 2 = %s" % (options.peaks2),\
            ""
            ))
    else:
        options.argtxt = "\n".join((
            "# ARGUMENTS LIST:",\
            "# name = %s" % (options.name),\
#            "# format = %s" % (options.format),\
            "# ChIP-seq file 1 = %s" % (options.t1bdg),\
            "# control file 1 = %s" % (options.c1bdg),\
            "# ChIP-seq file 2 = %s" % (options.t2bdg),\
            "# control file 2 = %s" % (options.c2bdg),\
            ""
            ))
         
        if options.peaks_pvalue:
            options.argtxt +=  "# treat/control -log10(pvalue) cutoff = %.2e\n" % (options.peaks_log_pvalue)
            options.argtxt +=  "# treat/control -log10(qvalue) will not be calculated and reported as -1 in the final output.\n"
        else:
            options.argtxt +=  "# treat/control -log10(qvalue) cutoff = %.2e\n" % (options.peaks_log_qvalue)
        
    if options.diff_pvalue:
        options.argtxt +=  "# differential pvalue cutoff = %.2e\n" % (options.log_pvalue)
        options.argtxt +=  "# differential qvalue will not be calculated and reported as -1 in the final output.\n"
    else:
        options.argtxt +=  "# differential qvalue cutoff = %.2e\n" % (options.log_qvalue)

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    return options

def opt_validate_filterdup ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logging.error("Error when interpreting --gsize option: %s" % options.gsize)
            logging.error("Available shortcuts of effective genome sizes are %s" % ",".join(list(efgsize.keys())))
            sys.exit(1)

    # format

    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser        
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # duplicate reads
    if options.keepduplicates != "auto" and options.keepduplicates != "all":
        if not options.keepduplicates.isdigit():
            logging.error("--keep-dup should be 'auto', 'all' or an integer!")
            sys.exit(1)

    # uppercase the format string 
    options.format = options.format.upper()

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    return options

def opt_validate_randsample ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # format

    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # uppercase the format string 
    options.format = options.format.upper()

    # percentage or number
    if options.percentage:
        if options.percentage > 100.0:
            logging.error("Percentage can't be bigger than 100.0. Please check your options and retry!")
            sys.exit(1)
    elif options.number:
        if options.number <= 0:
            logging.error("Number of tags can't be smaller than or equal to 0. Please check your options and retry!")
            sys.exit(1)

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    return options

def opt_validate_refinepeak ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # format

    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # uppercase the format string 
    options.format = options.format.upper()

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    return options

def opt_validate_predictd ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # gsize
    try:
        options.gsize = efgsize[options.gsize]
    except:
        try:
            options.gsize = float(options.gsize)
        except:
            logging.error("Error when interpreting --gsize option: %s" % options.gsize)
            logging.error("Available shortcuts of effective genome sizes are %s" % ",".join(list(efgsize.keys())))
            sys.exit(1)

    # format
    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
        options.nomodel = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
        options.nomodel = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "AUTO":
        options.parser = guess_parser
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # uppercase the format string 
    options.format = options.format.upper()

    # d_min is non-negative
    if options.d_min < 0:
        logging.error("Minimum fragment size shouldn't be negative!" % options.d_min)
        sys.exit(1)

    # upper and lower mfold
    options.lmfold = options.mfold[0]
    options.umfold = options.mfold[1]
    if options.lmfold > options.umfold:
        logging.error("Upper limit of mfold should be greater than lower limit!" % options.mfold)
        sys.exit(1)
    
    options.modelR  = os.path.join( options.outdir, options.rfile )

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    return options


def opt_validate_pileup ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # format

    options.gzip_flag = False           # if the input is gzip file
    
    options.format = options.format.upper()
    if options.format == "ELAND":
        options.parser = ELANDResultParser
    elif options.format == "BED":
        options.parser = BEDParser
    elif options.format == "ELANDMULTI":
        options.parser = ELANDMultiParser
    elif options.format == "ELANDEXPORT":
        options.parser = ELANDExportParser
    elif options.format == "SAM":
        options.parser = SAMParser
    elif options.format == "BAM":
        options.parser = BAMParser
        options.gzip_flag = True
    elif options.format == "BOWTIE":
        options.parser = BowtieParser
    elif options.format == "BAMPE":
        options.parser = BAMPEParser
        options.gzip_flag = True
    elif options.format == "BEDPE":
        options.parser = BEDPEParser
    else:
        logging.error("Format \"%s\" cannot be recognized!" % (options.format))
        sys.exit(1)
    
    # uppercase the format string 
    options.format = options.format.upper()

    # logging object
    logging.basicConfig(level=(4-options.verbose)*10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    # extsize
    if options.extsize <= 0 :
        logging.error("--extsize must > 0!")
        sys.exit(1)

    return options

def opt_validate_bdgcmp ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logging.basicConfig(level=20,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    # methods should be valid:

    for method in set(options.method):
        if method not in [ 'ppois', 'qpois', 'subtract', 'logFE', 'FE', 'logLR', 'slogLR', 'max' ]:
            logging.error( "Invalid method: %s" % method )
            sys.exit( 1 )

    # # of --ofile must == # of -m

    if options.ofile:
        if len(options.method) != len(options.ofile):
            logging.error("The number and the order of arguments for --ofile must be the same as for -m.")
            sys.exit(1)     

    return options


def opt_validate_cmbreps ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logging.basicConfig(level=20,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    # methods should be valid:


    if options.method not in [ 'fisher', 'max', 'mean']:
        logging.error( "Invalid method: %s" % options.method )
        sys.exit( 1 )

    if len( options.ifile ) < 2:
        logging.error("Combining replicates needs at least two replicates!")
        sys.exit( 1 )

    # # of -i must == # of -w

    # if not options.weights:
    #     options.weights = [ 1.0 ] * len( options.ifile )

    # if len( options.ifile ) != len( options.weights ):
    #     logging.error("Must provide same number of weights as number of input files.")
    #     sys.exit( 1 )

    # if options.method == "fisher" and len( options.ifile ) > 3:
    #     logging.error("NOT IMPLEMENTED! Can't combine more than 3 replicates using Fisher's method.")
    #     sys.exit( 1 )

    return options


def opt_validate_bdgopt ( options ):
    """Validate options from a OptParser object.

    Ret: Validated options object.
    """
    # logging object
    logging.basicConfig(level=20,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w"
                        )
    
    options.error   = logging.critical        # function alias
    options.warn    = logging.warning
    options.debug   = logging.debug
    options.info    = logging.info

    # methods should be valid:

    if options.method.lower() not in [ 'multiply', 'add', 'p2q', 'max', 'min']:
        logging.error( "Invalid method: %s" % options.method )
        sys.exit( 1 )

    if options.method.lower() in [ 'multiply', 'add' ] and not options.extraparam:
        logging.error( "Need EXTRAPARAM for method multiply or add!")
        sys.exit( 1 )

    return options

