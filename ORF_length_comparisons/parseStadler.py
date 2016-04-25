#150606 Mary K. Thompson
# parseStadler.py will 1) parse the DESeq files and 2) parse the gff files from wormbase
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
    '''Parse the gtf file obtained from wormbase'''

    f= open(gffFile, 'r')
    annd={}
    for line in f:
        if line.startswith('#'):
            continue
        elif line.startswith('###'):
            break
        fields= line.strip('\n').split('\t')
        chrom= fields[0]
        what= fields[1]
        region= fields[2]
        start= int(fields[3])
        end= int(fields[4])
        strand= fields[6]
        subattributes= fields[8].split(';')
        #changing this part to fit the UCSC format
        keys= [i.split('=')[0] for i in subattributes]
        vals=[]
        for i in subattributes:
            try:
                vals.append(i.split('=')[1]) #remove quotes
            except:
                vals.append('')
        attd= dict(zip(keys, vals))
    
        if (what=='Coding_transcript') and (region=='CDS'): #there are some things labeled CDS that don't have an associated protein for some reason
            txt= attd['ID'].split('CDS:')[1]
            protein= attd['wormpep'].split('CE:')[1]
           
            if protein not in annd: #could be there already if it's a second exon
                annd[protein]={}
                
            if txt not in annd[protein]:
                annd[protein][txt]={}
                annd[protein][txt]['exons']= [[start, end]]
                
            else:
                annd[protein][txt]['exons'].append([start, end])
    
    #get the length of the longest CDS for the protein
    newd={}
    for gene in annd:
        newd[gene]={}
        txts=[]
        for txt in annd[gene]:
            exonbounds= annd[gene][txt]['exons']
            exonbounds.sort(key=operator.itemgetter(0))
            
            exon_starts=[]
            exon_ends=[]
            
            for i in exonbounds:
                exon_starts.append(i[0])
                exon_ends.append(i[1])
            
            txts.append((exon_starts, exon_ends))
        
        longest=chooseLongest(txts)
        newd[gene]=longest
        
    return newd

def chooseLongest(txts):
    
    lens=[]
    for i in txts:
    
        ORFLen=0
        for j in range(0, len(i[0])):
            ORFLen+= i[1][j]- i[0][j]+1
        
        lens.append(ORFLen)
    
    longest= max(lens)
    return longest


def getTEs(FPfile, totFile, readco):
    '''
    Read the DESeq files and calculate TE as scaled FP reads/ scaled total reads. Make it so that both basemeans have >= readco
    Here A= before feeding, diapause; B= after feeding
    '''
    
    d={}
    f= open(FPfile, 'r')
    f.readline() #remove header
    for line in f:
        fields=line.strip('\n').split('\t')
        
        gene= fields[0]
        baseMeanA= float(fields[2])
        baseMeanB= float(fields[3])
        #expression_level= float(fields[2])
        #read_count= float(fields[3])
        
        d[gene]={}
        d[gene]['FP_A']= baseMeanA
        d[gene]['FP_B']= baseMeanB
    f.close()
    
    g= open(totFile, 'r')
    g.readline()
    for line in g:
        fields=line.strip('\n').split('\t')
        
        gene= fields[0]
        baseMeanA= float(fields[2])
        baseMeanB= float(fields[3])
        #expression_level= float(fields[2])
        #read_count= float(fields[3])
        if gene not in d:
            d[gene]={}
        d[gene]['tot_A']= baseMeanA
        d[gene]['tot_B']= baseMeanB
        
        #read_count could allow filtering later if needed
    g.close()
    
    print 'len d', len(d)
    #get TEs of all the genes:
    #this would be a good place to institute a read cutoff
    for gene in d:
        
        try:
            if (d[gene]['FP_A']+d[gene]['tot_A']>readco):
                TE_A= d[gene]['FP_A']/d[gene]['tot_A']
                d[gene]['TE_A']=math.log(TE_A, 10)
        except:
            pass
        
        try:
            if (d[gene]['FP_B']+d[gene]['tot_B']>readco):
                TE_B= d[gene]['FP_B']/d[gene]['tot_B']
                d[gene]['TE_B']=math.log(TE_B, 10)
        except:
            pass
    
    return d

def correlate(lenD, TEs, outname):
    '''report spearman r between ORF length and transcript length'''
    xA=[]
    yA=[]
    newdA={}
    for gene in lenD:
        try:
            x1= lenD[gene]
            y1= TEs[gene]['TE_A']
            xA.append(x1)
            yA.append(y1)
            newdA[gene]={}
            newdA[gene]['TE']=y1
            newdA[gene]['ORFLen']=x1
        except:
            continue
    
    print 'A'
    print 'lenx', len(xA)
    spear= stats.spearmanr(xA, yA)
    print 'spear', spear
    
    xB=[]
    yB=[]
    newdB={}
    for gene in lenD:
        try:
            x1= lenD[gene]
            y1= TEs[gene]['TE_B']
            xB.append(x1)
            yB.append(y1)
            newdB[gene]={}
            newdB[gene]['TE']=y1
            newdB[gene]['ORFLen']=x1
        except:
            continue
    
    print 'B'
    print 'lenx', len(xB)
    spear= stats.spearmanr(xB, yB)
    print 'spear', spear
            
   
    f=open('%s_A.p' % outname, 'w')
    pickle.dump(newdA, f)
    f.close()
    
    f=open('%s_B.p' % outname, 'w')
    pickle.dump(newdB, f)
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
        lenD= parseGtf(args['annotations'])
        f= open('%s_len.p' % args['outname'], 'w')
        pickle.dump(lenD, f)
        f.close()
        sys.exit() #quit b/c assume don't want to run the rest for now
    else:
        f= open(args['annotations'], 'r')
        lenD= pickle.load(f)
        f.close()
        
    readco=128 #must have at least this many scaled reads to calculate TE
    
    #check that this looks good:
    div3=0
    notdiv3=0
    allLens=[]
    for gene in lenD:
        length=lenD[gene]
        allLens.append(length)
        if length%3==0:
            div3+=1
        else:
            notdiv3+=1
    
    print 'div3', div3
    print 'notdiv3', notdiv3
    print 'meanlen', np.mean(allLens)
    
    TEs= getTEs(args['FP_expFile'], args['tot_expFile'], readco)
    correlate(lenD, TEs, args['outname']) #get TE vs. length relationship and save as output file in format d[gene][length]=length; d[gene][TE]=TE
    
if __name__ == '__main__':
    main(sys.argv[1:])
    