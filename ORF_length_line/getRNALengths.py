#150507 Mary K. Thompson
#purpose: put yeast ORF vs. transcript lengths into a convenient format for further manipulation

import sys
import operator
import numpy as np
import cPickle as pickle

pelechanoFile= '/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/C6.24/nature12121-s2/Supplementary Data 3.txt'
pelechanoIsoforms='/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/C6.24/nature12121-s2/S2_tcd_mTIFAnno.txt'
yeast_gff='/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/C6.24/SacCer3_annotations/20110902_sacCer3.gff'

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

def parsePelechano(pelechanoFile):
    
    pelechanoD={}
    f= open(pelechanoFile, 'rU')
    header= f.readline().strip('\n').split('\t')
    for line in f:
        fields=line.strip('\n').split('\t')
        geneID= fields[2]
        pelechanoD[geneID]={}
        med5= float(fields[5])
        med3= float(fields[6])
        pelechanoD[geneID]['med5']= med5
        pelechanoD[geneID]['med3']= med3
    
    f.close()
    return pelechanoD

def parsePelechanoIsoforms(pelechanoIsoforms):
    '''
    parse the major transcript isoforms file because they have matched 5' and 3' boundaries so that I can use this to get the
    average transcript length for each gene, weight the numbers by including only the ypd present transcripts
    '''
    f=open(pelechanoIsoforms, 'rU')
    header= f.readline().strip('\n').split(' ')
    isoformD={}
    
    for line in f:
        fields= line.strip('\n').split(' ')
        start= int(fields[2])
        stop= int(fields[3])
        length= abs(start-stop)+1 #I'm assuming that this is endpoint inclusive, but would be good to check
        ypd= int(fields[4])
        gal= int(fields[5])
        description= ' '.join(fields[6:])
        
        if 'Covering one intact ORF' in description: #only add ones that completely cover an ORF
            ORFName= fields[-1]
            if ypd >0:
                if ORFName not in isoformD:
                    isoformD[ORFName]={}
                isoformD[ORFName][length]= ypd
    f.close()
    #now get average and median lengths in ypd using this data
    avIsoD= averageIsoforms(isoformD)
    return avIsoD

def averageIsoforms(isoformD):
    '''
    return average and median transcript isoform length from the Pelechano data
    '''
    isoD={}
    for gene in isoformD:
        isoD[gene]={}
        x=[]
        for i in isoformD[gene]:
            x.extend([i]*isoformD[gene][i]) #adds the number the appropriate number of times
        
        medLen= np.median(x)
        avLen= np.mean(x)
        isoD[gene]['median']= medLen
        isoD[gene]['average']= avLen
    
    return isoD
    
def getTxts(txtboundaryD, isoformD, ORFLenD):
    '''
    combine the lengths of the ORFs with their median 5' and 3' lengths from Pelechano. Note this is not really a perfect analysis
    because these median values are not matched (i.e. maybe the longer 5' UTRs generally go with the shorter 3' UTRs, etc.).
    '''
    
    txtsD={}
    
    for gene in ORFLenD:
        txtsD[gene]={}
        if gene in txtboundaryD:
            totalLen= ORFLenD[gene]['ORFLen']+txtboundaryD[gene]['med5']+txtboundaryD[gene]['med3']
    
            txtsD[gene]['sepLen']= totalLen
    
        if gene in isoformD:
            txtsD[gene]['median']= isoformD[gene]['median']
            txtsD[gene]['average']= isoformD[gene]['average']
            
        txtsD[gene]['ORFLen']= ORFLenD[gene]['ORFLen']
        
    return txtsD

def main():
    
    yeastOutname= sys.argv[1]
    
    #yeast annotations:
    txtboundaryD= parsePelechano(pelechanoFile)
    ORFLenD= getORFLens(yeast_gff)
    isoformD= parsePelechanoIsoforms(pelechanoIsoforms)
    yeasttxtsD= getTxts(txtboundaryD, isoformD, ORFLenD)
    
    f= open('%s.p' % yeastOutname, 'w')
    pickle.dump(yeasttxtsD, f)
    f.close()
    
main()
