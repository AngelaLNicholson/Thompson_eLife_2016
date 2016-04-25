#150606 Mary K. Thompson
#parseGuo.py will 1) parse the rpkm files and 2) parse the gff files and return transcript lengths
#since Guo actually reports the transcript it's mapped to, I will use those instead of the longest
import sys
import argparse
import operator
import numpy as np
import math
import scipy.stats as stats
import cPickle as pickle

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

def parseGtf(gffFile):
    '''Parse the gtf file obtained from UCSC'''
    f= open(gffFile, 'r')
    annd={}
    for line in f:
        if line.startswith('#'):
            continue
        elif line.startswith('###'):
            break
        fields= line.strip('\n').split('\t')
        chrom= fields[0]
        region= fields[2]
        start= int(fields[3])
        end= int(fields[4])
        strand= fields[6]
        subattributes= fields[8].split('; ')
        #changing this part to fit the UCSC format
        keys= [i.split(' ')[0] for i in subattributes]
        vals=[]
        for i in subattributes:
            try:
                vals.append(i.split(' ')[1].strip('"')) #remove quotes
            except:
                vals.append('')
        attd= dict(zip(keys, vals))
    
        if region=='CDS':
            ID= attd['transcript_id']
            if ID not in annd: #could be there already if it's a second exon
                annd[ID]={}
                annd[ID]['chr']= chrom
                annd[ID]['strand']= strand
                   
                annd[ID]['exons']= [[start, end]]
                
            else:
                annd[ID]['exons'].append([start, end])
    
    #sort the exons by boundary, assume that there can only be one leftmost exon and one rightmost. the rightmost will also be the end of the gene
    for gene in annd:
        exonbounds= annd[gene]['exons']
        exonbounds.sort(key=operator.itemgetter(0))
        
        annd[gene]['start']= exonbounds[0][0]
        annd[gene]['end']= exonbounds[-1][1]
        
        exon_starts=[]
        exon_ends=[]
        
        for i in annd[gene]['exons']:
            exon_starts.append(i[0])
            exon_ends.append(i[1])
            
        annd[gene]['exon_starts']= exon_starts
        annd[gene]['exon_ends']= exon_ends
        
    return annd


def getORFLens(gffFile):
    annd= parseGtf(gffFile)
    
    lenD={}
    
    for gene in annd:
        lenD[gene]={}
        ORFLen=0
        for i in range(0, len(annd[gene]['exons'])):
            ORFLen+= annd[gene]['exons'][i][1]- annd[gene]['exons'][i][0]+1
        
        lenD[gene]['ORFLen']= ORFLen
    
    return lenD

def getTEs(FPfile, totFile, readco):
    
    
    d={}
    f= open(FPfile, 'r')
    f.readline() #remove header
    for line in f:
        fields=line.strip('\n').split('\t')
        
        gene= fields[0]
        expression_level= float(fields[2])
        read_count= float(fields[3])
        
        d[gene]={}
        d[gene]['FP_rpkm']= expression_level
        d[gene]['FP_counts']= read_count
    f.close()
    
    g= open(totFile, 'r')
    g.readline()
    for line in g:
        fields=line.strip('\n').split('\t')
        
        gene= fields[0]
        expression_level= float(fields[2])
        read_count= float(fields[3])
        if not gene in d:
            d[gene]={}
        d[gene]['tot_rpkm']= expression_level
        d[gene]['tot_counts']= read_count
        #read_count will allow filtering later
    g.close()
    
    #get TEs of all the genes:
    #this would be a good place to institute a read cutoff
    for gene in d:
        try:
            if (d[gene]['FP_counts'] + d[gene]['tot_counts'])>= readco:
                TE= d[gene]['FP_rpkm']/d[gene]['tot_rpkm']
                d[gene]['TE']=math.log(TE, 10)
            else:
                continue
        except:
            continue
        
    return d

def correlate(lenD, TEs, outname):
    '''report spearman r between ORF length and transcript length'''
    x=[]
    y=[]
    newd={}
    for gene in lenD:
        try:
            x1= lenD[gene]['ORFLen']
            y1= TEs[gene]['TE']
            x.append(x1)
            y.append(y1)
            newd[gene]={}
            newd[gene]['TE']=y1
            newd[gene]['ORFLen']=x1
        except:
            continue
    
    print 'lenx', len(x)
    spear= stats.spearmanr(x, y)
    print 'spear', spear
    
    f=open('%s.p' % outname, 'w')
    pickle.dump(newd, f)
    f.close()    
    
def main(argList):
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('annotations')
    parser.add_argument('FP_expFile')
    parser.add_argument('tot_expFile')
    parser.add_argument('outname')
    parser.add_argument('--writeLenP', action='store_true')
   
    ar=parser.parse_args(args=argList)
    args= vars(ar)
    
    if args['writeLenP']==True: #set this flag if you want to write the lens.p file, otherwise provide the lens.p file in place of the gtf file in the arglist
        lenD= getORFLens(args['annotations'])
        f= open('%s_len.p' % args['outname'], 'w')
        pickle.dump(lenD, f)
        f.close()
        sys.exit() #quit b/c assume don't want to run the rest for now
    else:
        f= open(args['annotations'], 'r')
        lenD= pickle.load(f)
        f.close()
    
    #check that this looks good:
    div3=0
    notdiv3=0
    allLens=[]
    for gene in lenD:
        length=lenD[gene]['ORFLen']
        allLens.append(length)
        if length%3==0:
            div3+=1
        else:
            notdiv3+=1
    
    print 'div3', div3
    print 'notdiv3', notdiv3
    print 'meanlen', np.mean(allLens)
    
    
    readco=128
    
    TEs= getTEs(args['FP_expFile'], args['tot_expFile'], readco)
   
    correlate(lenD, TEs, args['outname']) #get TE vs. length relationship and save as output file in format d[gene][length]=length; d[gene][TE]=TE
    
    
if __name__ == '__main__':
    main(sys.argv[1:])
    