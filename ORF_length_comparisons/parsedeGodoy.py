#150609 Mary K. Thompson
#purpose: parse mass spec data from de Godoy et al., 2008 and compare with deltaTE in asc1-M1X

import sys
import csv
import scipy.stats as stats
import cPickle as pickle

def read(infile):
    
    f= open(infile, 'rU')
    reader= csv.reader(f)
    header= reader.next()
    
    lightI= header.index('Intensity L')
    heavyI= header.index('Intensity H')
    expI= header.index('Experiment') #this will either be 'Lys-C' or 'Trypsin'

    exps=set()
    genes=[]
    d={}
    for row in reader:
        try:
            sysname= row[0]
            lightval= float(row[lightI])
            heavyval= float(row[heavyI])
            exptype= row[expI]
            exps.add(exptype)
            genes.append(sysname)
            
            if exptype=='Trypsin':
                d[sysname]= heavyval
        
            else:
                d[sysname]= lightval
                 
        except:
            continue
    f.close()
   
    #check if there are replicates in the set of genes
    seen = set()
    uniq = []
    for x in genes:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    
    print 'len seen', len(seen)
    print 'len unique', len(uniq) #it looks like this is unique, only Lys-C or Trypsin is reported/protein
    
    return d
    
def correlate(quant, d):
    '''
    get corrleation between delta TE and msquant data
    '''
    x=[]
    y=[]
    
    for gene in d:
        if 'TE' in d[gene]:
            try:
                x1= d[gene]['TE']['val']
                y1= quant[gene]
                x.append(x1)
                y.append(y1)
            except:
                continue
    
    spearman= stats.spearmanr(x, y)
    print 'spearman', spearman
           
def main():
    infile, datapickle = sys.argv[1:]
    
    msquant=read(infile)
    
    f= open(datapickle, 'r')
    d= pickle.load(f)
    f.close()
    
    correlate(msquant, d)
    

main()