#141028 Mary K. Thompson
#Purpose: generate figure showing (mRNA attributes vs. M1X delta TEs), add readcount filtering to the WT TE part

import sys
import csv
import cPickle as pickle
import operator
import numpy as np
import pandas as pd
import math

def parseGff(gffFile):
    '''Parse the gff'''
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
        subattributes= fields[8].split(';')
        keys= [i.split('=')[0] for i in subattributes]
        vals=[]
        for i in subattributes:
            try:
                vals.append(i.split('=')[1])
            except:
                vals.append('')
        attd= dict(zip(keys, vals))
        #print attd
    
        if region=='CDS':
            ID= attd['Parent'].split('_')[0] #get rid of the _mRNA part
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
    annd= parseGff(gffFile)
    
    lenD={}
    
    for gene in annd:
        lenD[gene]={}
        ORFLen=0
        for i in range(0, len(annd[gene]['exons'])):
            ORFLen+= annd[gene]['exons'][i][1]- annd[gene]['exons'][i][0]+1
        
        lenD[gene]['ORFLen']= ORFLen
    
    return lenD

def parseCsv(csvFile):
    f= open(csvFile, 'rU')
    reader= csv.reader(f)
    header= next(reader)
    
    attD={}
    for row in reader:
        attD[row[0]]={}
        for i in range(1, len(header)):
            attD[row[0]][header[i]]= row[i]
    f.close()
    return attD

def parsePelechano(infile):
    
    pUTRlens={}
    
    f= open(infile, 'r')
    header= f.readline().strip('\n').split('\t')
    for line in f:
        fields= line.strip('\n').split('\t')
        genename= fields[2]
        pUTRlens[genename]={}
        pUTRlens[genename]['median5']= int(round(float(fields[5])))
        pUTRlens[genename]['median3']= int(round(float(fields[6])))
    return pUTRlens
    
def getMFEs(MFE5p, MFE3p):
    '''
    given the files containing the calculated 5p and 3p MFEs, assimilate into a dictionary
    '''
    energyD={}
    
    f= open(MFE5p, 'r')
    g= open(MFE3p, 'r')
    
    for line in f:
        fields= line.strip('\n').split('\t')
        gene= fields[0].split('_')[0]
        fMFE= float(fields[1])
        
        energyD[gene]={}
        energyD[gene]['5pMFE']= fMFE
    f.close()
    
    for line in g:
        fields= line.strip('\n').split('\t')
        ID= fields[0].split('_')[0]
        tMFE= float(fields[1])
        
        if ID not in energyD:
            energyD[ID]={}
            
        energyD[ID]['3pMFE']= tMFE
        
    g.close()
    return energyD
    
def main():
    attfile, pelechanoUTRFile, MFE5p, MFE3p, gfffile, datapickle, unfilteredpickle, rpkmpickle, outprefix= sys.argv[1:]
    
    bigD={}
    
    #construct large dictionary {gene:'deltaTE':___, '5pUTR':__, '3pUTR':___, etc. }
    f= open(datapickle, 'r')
    d= pickle.load(f)
    f.close()
    
    g= open(unfilteredpickle, 'r')
    unf= pickle.load(g)
    g.close()
    
    #populate the dictionary with all of the genes that we have deltaTE values for
    for gene in d:
        if 'TE' in d[gene]:
            bigD[gene]={}
            bigD[gene]['deltaTEM1X']= d[gene]['TE']['val']
    
   
    #the rpkm calculations take the mitochondrial genome values into account in the total whereas the scaled counts do not.
    #there are some genes present in the rpkm dict but not in the filtered values b/c the WTTE=0
    
    g= open(rpkmpickle, 'r')
    rd= pickle.load(g)
    g.close()
    for gene in rd:
        TEreads=sum([rd[gene]['WT_T1']['counts'], rd[gene]['WT_T2']['counts'], rd[gene]['WT_Fp1']['counts'], rd[gene]['WT_Fp2']['counts']])
        if TEreads<128:
            continue
        else:
            wtrpkm= math.log(np.mean([rd[gene]['WT_T1']['rpkm'], rd[gene]['WT_T2']['rpkm']]), 10)
            try:
                wtTE= math.log(unf[gene]['TE']['baseMeanA'], 10)
            except:
                wtTE='NaN'
            if gene not in bigD:
                bigD[gene]={}
            bigD[gene]['WT_TE']= wtTE
            bigD[gene]['WT_rpkm']= wtrpkm
                
    #parse the other contributing datasets:
    attD= parseCsv(attfile)
    pelechanoUTRs= parsePelechano(pelechanoUTRFile)
    ORFLenD= getORFLens(gfffile)
    MFEs= getMFEs(MFE5p, MFE3p)
    
    #add the attributes that we're interested in from each dataset:    
    atts2add={'attD':['tAI', 'pALe'], 'pelechanoUTRs':['median3', 'median5'], 'ORFLenD':['ORFLen'], 'MFEs':['5pMFE', '3pMFE']}
    stringNames={'attD':attD, 'pelechanoUTRs':pelechanoUTRs, 'ORFLenD':ORFLenD, 'MFEs':MFEs}
    
    #for gene in bigD:
    for gene in ORFLenD: #now it will add attributes for everything in the annotations
        if gene not in bigD:
            bigD[gene]={}
        for s in atts2add:
            for key in atts2add[s]:
                try:
                    bigD[gene][key]= stringNames[s][gene][key]
                except:
                    bigD[gene][key]='NaN'
    
    
    bigDataFrame= pd.DataFrame.from_dict(bigD, orient='index')
    bigDataFrame.to_csv(path_or_buf='%s.csv' % outprefix)
       
main()