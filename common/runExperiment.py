#140616 Mary K. Thompson
#purpose: parse the input file, run the analysesRunner and graphPlotters scripts

import sys
import os
import argparse
from collections import defaultdict
import cPickle as pickle
import math
import graphPlotters as graphPlotters
import analysesRunner as analysesRunner

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
        self.getFilePaths()
        
    def parseParams(self, paramfile):
        
        '''
        will list INDIR and OUTDIR separately
        allow a number after IND1R and OUTDIR, i.e. INDIR1, so that we can just use the os module to create the full path names
        for PLOTS, will save the plotting options as a list so that it can just be given as an arglist to the graph_plotters module
        '''
        block='' #keep track of what block we're in
        num='' #keep track of any numbers associated with the block
        
        self.directories={}
        
        self.directories['indirs']=defaultdict() #use num->list of directories or files
        self.directories['outdir']={} #use num->list of directories or files        
        self.files=defaultdict(list)
        
        self.analyses=[] #list to store the Analysis objects that will be used to run specific ones later
        self.plots=[]
        #self.keys={} #dictionary to store the keys to compare, i.e. keys={'WT':'baseMeanA', 'M1X':'baseMeanB'}
        self.outprefix=''
        
        f= open(paramfile, 'r')
        #d= defaultdict()
        
        for line in f:
            line= line.strip('\n')
           
            if (len(line)==0) or (line.startswith('#')):
                continue
            
            elif line.startswith('INDIR'): #ignore whitespace on either side of the equals sign
                fields= stripLine(line, '=')
                dirnum= int(fields[0].split('INDIR')[1])
                indir= fields[1]                  
                self.directories['indirs'][dirnum]=indir #if I do this, then we'll negate the purpose of having an indir name
            
            elif line.startswith('OUTDIR'):
                fields= stripLine(line, '=')
                outdir=fields[1]
                self.directories['OUTDIR']=outdir #if I do this, then we'll negate the purpose of having an indir name
                #note that I'm assuming for now that all output will go to one directory, therefore they are all in the '' key

            elif line.startswith('OUTPREFIX'):
                fields= stripLine(line, '=')
                outprefix= fields[1]
                self.outprefix= outprefix
                
            elif line.startswith('FILES'):
                block='file'
                num= int(line.split('FILES')[1])
                
            elif line.startswith('ANALYSES'):
                block= 'analyses'
                
            elif line.startswith('PLOTS'):
                block= 'plots'
           
            elif block== 'file':                
                filename= line
                self.files[num].append(filename)
            
            elif block=='analyses':
                fields= stripLine(line, ' ')
                analysisargs= self.replaceDirs(fields)
                self.analyses.append(analysisargs)

            elif block=='plots':
                fields= stripLine(line, ' ')
                plotargs= self.replaceDirs(fields)
                self.plots.append(plotargs)
                #in order to have lines that can have names but no keywords (makes easiest to edit param file this way), you need to have truth blocks and keep track of what block you're in
        f.close()
    
    def replaceDirs(self, fields):
        '''Replace any variables-- i.e. directory names like {OUTDIR} with the actual name'''
        for i in range(0, len(fields)):
            if '{' in fields[i]:
                parts= fields[i].split('}')
                base= parts[1] #get the thing to the right of the brackets, assuming that we're only using this for directory names
                toreplace= parts[0][1:] #get the name free of the {} brackets
                fulldir=self.directories[toreplace]
                newdir=''.join([fulldir, base])
                fields[i]= newdir
                
        return fields
        
    def getFilePaths(self):
        '''
        Replace files with fullnames, i.e. including the indir
        '''
        newfiles=[]
        thesedirs= self.files.keys()
        thesedirs.sort()
        for i in thesedirs:
        #for num in self.files:
            for shortname in self.files[i]:
                newname= os.path.join(self.directories['indirs'][i], shortname)
                newfiles.append(newname)
        
        self.files= newfiles #reassign to just contain file names
                
def main(argList):
        
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('paramfile', help='paramfile describing which analyses and plots to run')
        
    ar= parser.parse_args(args=argList)
    args= vars(ar)
    
    paramfile= args['paramfile']
    params= Params(paramfile)
    if not os.path.exists(params.directories['OUTDIR']):
        os.makedirs(params.directories['OUTDIR'])
    
    for analysis in params.analyses:
        analysis.extend(['-infiles'])
        analysis.extend(params.files)
        print 'analysis', analysis
        #RIGHT NOW WE CAN MAKE A LIST OF OBJECTS AND PICKLE IT OR WE CAN JUST MAKE A SINGLE OBJECT AND PICKLE IT. LIST WOULD BE USED FOR SOMETHING LIKE BOX PLOT WITH DIFFERENT MUTANTS AS THE SERIES
        if '-list' in analysis: #if we want it to be a sequential list of objects now going to have to specify with the --list option
            obj_list=[]
            for i in range(0, len(params.files)):
                analysis.extend(['-filenum'])
                analysis.extend([i])
                obj, aoutname= analysesRunner.main(analysis)
                obj_list.append(obj)
            obj= obj_list #assign so that the list will get pickled
        
        else:
            obj, aoutname= analysesRunner.main(analysis)

        #pickle the object to use for graphing later
        if params.outprefix!='':    
            basename='%s_%s' % (params.outprefix, aoutname)
        else:
            basename= aoutname
        
        outname= '%s.p' % os.path.join(params.directories['OUTDIR'],basename)
        
        print 'pickling obj', obj   
        f=open(outname, 'w')
        pickle.dump(obj, f)
        f.close()
        
    for plot in params.plots:
        graphPlotters.main(plot)
    
if __name__ == '__main__':
    main(sys.argv[1:])
