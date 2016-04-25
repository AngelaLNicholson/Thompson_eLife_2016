#140611 Mary K. Thompson
#purpose: wrapper script to run the scripts required to go from the CDSprints.p files to the experiment objects used for downstream analyses and plots

#import python utilities:
import sys
import os
from collections import defaultdict
import cPickle as pickle
import math
import argparse
import subprocess
import combineAllDEdata
import combineAllDEdata_noReplicates

#import scripts needed for DESeq-related processing:
import CDSprints2DESeq
import filterByCount
import makeGEOfiles

def stripLine(line, splitter):
    '''
    returns the components of a line with all whitespaces removed and split on splitter
    '''
    l=[]
    fields= line.split(splitter)
    for i in fields:
        l.append(i.strip())
    
    return l

class Params(object):
    
    params={}
    def __init__(self, file):
        self.parseParams(file)
        
    def parseParams(self, paramfile):
        
        '''
        will list INDIR and OUTDIR separately
        allow a number after IND1R and OUTDIR, i.e. INDIR1, so that we can just use the os module to create the full path names
        for PLOTS, will save the plotting options as a list so that it can just be given as an arglist to the graph_plotters module
        '''
        block='' #keep track of what block we're in
        num='' #keep track of any numbers associated with the block
        
        f= open(paramfile, 'r')
        
        self.libraries={} #maps [libname]->full path
        self.replicates={} #maps [strain_name]=[FP_r1, FP_r2, T_r1, T_r2]
        self.bmcnames={} #maps libname-> bmc .fastq.gz name
        self.lib_order=[] #this will be the for each library, i.e. WT_Fp1, WT_Fp2, ..
        self.extra_order=[]
        self.strain_order=[] #this will be WT, D109Y, etc.
        self.reporder={} #store the strain names matched with the WT control they belong to, i.e., {0:[WT, D109Y, M1X, asc1Null], 1:[WT_pA, DE_pA]}
        self.WT=[] #store the name of the matched WT control
        
        for line in f:
            line= line.strip('\n')
           
            if (len(line)==0) or (line.startswith('#')):
                continue
            
            elif line.startswith('summary_outname'): #ignore whitespace on either side of the equals sign
                fields= stripLine(line, '=')
                summary_outname= fields[1]
                self.bigoutname= summary_outname
            
            elif line.startswith('experiment_outdir'):
                fields= stripLine(line, '=')
                outdir=fields[1]
                self.outdir=outdir #if I do this, then we'll negate the purpose of having an indir name
                #note that I'm assuming for now that all output will go to one directory, therefore they are all in the '' key

            elif line.startswith('codons_to_exclude'):
                fields= stripLine(line, '=')
                codons= fields[1]
                self.codons= int(codons)
            
            elif line.startswith('read_co'):
                fields= stripLine(line, '=')
                co= fields[1]
                self.co= int(co)
            
            elif line.startswith('min_map_len'):
                fields= stripLine(line, '=')
                minmaplen= fields[1]
                self.minmaplen= int(minmaplen)
                
            elif line.startswith('WT_condition'):
                fields=stripLine(line, '=')
                condition= fields[1]
                self.WT.append(condition)
                self.reporder[len(self.WT)-1]=[]
            
            elif line.startswith('GEO_OUTDIR'):
                fields=stripLine(line, '=')
                self.geo_outdir= fields[1]
            
            elif line.startswith('BMC_dir'):
                fields=stripLine(line, '=')
                self.bmcdir= fields[1]
                
            elif line.startswith('library_locations'):
                block='library'
                
            elif line.startswith('extra_libraries'):
                block='extra'
            
            elif line.startswith('replicate_matching'):
                block= 'replicates'
            
            elif line.startswith('BMC_Names'):
                block= 'bmc'
                
            elif block== 'library':                
                fields= stripLine(line, '=')
                name= fields[0]
                self.lib_order.append(name)
                path= fields[1]
                self.libraries[name]= path
                
            elif block=='extra':
                fields= stripLine(line, '=')
                name= fields[0]
                self.extra_order.append(name)
                
            elif block=='replicates':
                fields= stripLine(line, '=')
                #split the names of the libraries individually
                libs= fields[1].split(' ')
                self.replicates[libs[0]]= libs[1:]
                self.reporder[len(self.WT)-1].append(libs[0])
                self.strain_order.append(libs[0])
            
            elif block=='bmc':
                fields= stripLine(line, ' ')
                bmcbase= fields[0]
                libname= fields[1]
                self.bmcnames[libname]=bmcbase
        
        f.close()
        
def main(argList):
        
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('paramfile', help='parameter file describing file locations and analysis preferences')
        
    ar= parser.parse_args(args=argList)
    args= vars(ar)
    
    print 'parsing Params'
    #get parameters and name-matching
    paramfile= args['paramfile']
    params= Params(paramfile)
    
    print 'running CDSprints2DESeq'
    #get rpkms and count data    
    #also create the Rscript file needed for DESeq
    CDSprints2DESeq.main([params])
    
    print 'running R'
    #run DESeq via R
    Rscriptname= os.path.join(params.outdir, params.bigoutname)+'.R'
    cmd= '%s %s' % ('Rscript', Rscriptname)
    subprocess.check_call(cmd, shell=True)
    
    print 'combining all data into experiment dictionaries'
    #generate the _allData.p experiment files for each experiment
    combineAllDEdata.main([params])
    
    print 'writing GEO files, including md5sums'
    #write GEO output files, get md5sums as well
    makeGEOfiles.main([params])
    
    print 'filtering by read count to write filtered files'
    #create read count filtered files using 128 total counts as total of all libraries cutoff to create the _filtered.p experiment files
    filterByCount.main([params])

if __name__ == '__main__':
    main(sys.argv[1:])