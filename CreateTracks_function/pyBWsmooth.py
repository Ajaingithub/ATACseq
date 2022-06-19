#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University
# This will make a bigWig of ATAC data

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
import subprocess
from optparse import OptionParser
import pysam
from multiprocessing import Pool
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile
# http://bcbio.wordpress.com/2009/04/29/finding-and-displaying-short-reads-clustered-in-the-genome
# http://nullege.com/codes/search/bx.intervals.Intersecter?fulldoc=1  

#### OPTIONS ####
# read options from command line
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-p",default=20,type='int', help="<Threads> Accepts an integer")
opts.add_option("-a", help="<bw> Accepts a bigwig file")
opts.add_option("-g", help="<Genome Size file>")
opts.add_option("-w", default=150,type='int', help="<Int> window size")
opts.add_option("-s", default=20,type='int', help="<Int> step size (span)")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### DEFINE FUNCTIONS ##### 

##### INPUTS AND OUTPUTS #####
# open bigwig
bw = BigWigFile(open(options.a))

# get gSize file
gSizes = np.loadtxt(options.g,'str')
chunkSize = 1000000
padLen = 5000

# open out file
outName = options.a.split('.dgf.bw')[0]+'out.smooth.bed'
try: os.remove(outName)
except OSError: pass
outF = file(outName, 'a')

# smoothParems
wSize = options.w
wSmooth = np.ones(wSize)
step = options.s

#### SCRIPT #####
# split genome into chunks
for i in range(0,len(gSizes)):
    # break chrs into pieces
    chrN = gSizes[i][0]
    chrLen = int(gSizes[i][1])
    sVals = np.arange(1,int(gSizes[i][1]),chunkSize)
    
    # read in bigWig
    for j in range(0,len(sVals)):
        # get data, pass if not available
        signal = bw.get_as_array(chrN,sVals[j],sVals[j]+chunkSize+padLen)
        try: signal.any()
        except: continue
        
        # smooth data
        print chrN, sVals[j]
        signal[np.isnan(signal)] = 0
        convM = np.convolve(signal,wSmooth,'same')
        
        # save data
        sList = np.arange(sVals[j],sVals[j]+chunkSize+padLen,step)
        eList = sList+step
        chrList = np.array([chrN]*len(sList))
        meanSig = convM[range(step/2,chunkSize+padLen,step)]
        
        # save out
        idx1 = meanSig>0; idx2 = eList<chrLen; idx = idx1*idx2
        idx[chunkSize/step:] = False
        pData = np.c_[chrList[idx],np.array(sList[idx],dtype=str),eList[idx],meanSig[idx]]
        np.savetxt(outF, pData, fmt='%s',delimiter='\t', newline='\n')

# close outF
outF.close()

# convert to bw
print 'Converting to bigwig...'
out = options.a.split('.dgf.bw')[0]+'.s'+str(step)+'.w'+str(wSize)+'sw.bw'
tool='/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/CreateTracks/bedGraphToBigWig'
os.system(tool+' '+outName+' '+options.g+' '+out)

# sign out
print "Created by Jason Buenrostro."
print "Completed"
