#140317 Mary K. Thompson
#purpose: Perform analyses on experiment (_allData.p) files and store output for later plotting

import sys
import os
import matplotlib
import argparse
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import numpy as np
import scipy.stats as stats
import copy

from scipy.stats import gaussian_kde
from numpy.random import normal
from numpy import arange
from matplotlib.patches import Polygon
from collections import defaultdict
import cPickle as pickle
import math
import graphPlotters as graphplotters

infinity= float('inf')
taboo=set([infinity, -infinity, 'ND'])

def unpickle(infile):
    f= open(infile, 'r')
    d= pickle.load(f)
    f.close()
    return d

class histObject(object):
    
    def __init__(self, allCdfs, allHists, allEdges, pvals, legends, genecounts, allVals, mutants=[]):
        self.allCdfs= allCdfs
        self.allHists= allHists
        self.allEdges= allEdges
        self.pvals= pvals
        self.legends= legends
        self.mutants= mutants
        self.genecounts= genecounts
        self.allVals= allVals

class scatterObject(object):
    
    def __init__(self, allPlotData, spear_rvalues, spear_pvalues, pear_rvalues, pear_pvalues, legends=[], genecounts=[], order=[], mutants=[]):
        self.allPlotData= allPlotData
        self.spear_rvalues= spear_rvalues
        self.spear_pvalues= spear_pvalues
        self.pear_rvalues= pear_rvalues
        self.pear_pvalues= pear_pvalues
        self.legends= legends
        self.order=order
        self.mutants= mutants
        self.genecounts=genecounts
        self.flipAxes= False
        self.otype='Scattero'

class SingleScatterObject(object):
    '''
    If given a multiscatter object and coordinates, will create a single scatter object-- to use all the same code.
    '''
    
    def __init__(self, d, pair):
            
        index1= pair[0][0]
        index2= pair[0][1]
        y= d.allPlotData[index1]
        x= d.allPlotData[index2] #b/c of the way I did the counting, actually need to swap x and y like that

        datasets=[[x,y]]
        
        self.allPlotData= datasets
        self.legends= [d.legends[index1], d.legends[index2]]
        self.flipAxes= False
        self.otype='Histo'
            

def addGenes(geneSetFiles):
    
    if geneSetFiles==None:
        return []
    taboo=set(['>', ' ', '#']) #skip the line if starts with any of these
    
    sets=[]
    
    for i in geneSetFiles:
        geneset=set()
        f= open(i, 'r')
        for line in f:
            if line[0] in taboo:
                continue
            gene=line.strip('\n')
            geneset.add(gene)
        f.close()
        sets.append(geneset)
    return sets
    
def getVals(geneset='all'):
    '''
    return values for genes, default geneset is all genes for which data is available
    '''
    
    comparex= args['comparex']
    comparey= args['comparey']
    datasets= args['datasets']
    
    #figure out which values to collect for the x and yvalues:
    xdataset= int(comparex[0])-1
    xdatatype= comparex[1]
    xdatakey= comparex[2]
    
    #TWO OPTIONS SUPPORTED, GO THROUGH IN LOOP IF 'LIST' OR GET ONES SPECIFICALLY FROM THAT DATASET
    if geneset=='all':
        if args['list']!= None:
            geneset= set(datasets[0].keys())
            for j in range(1, len(datasets)):
                geneset=geneset.intersection(set(datasets[j].keys()))
        else:
            geneset= set(datasets[xdataset].keys())
    
    vals=[]
    
    for gene in geneset:
        if gene in datasets[xdataset]:
            if xdatatype in datasets[xdataset][gene]:
                if xdatakey in datasets[xdataset][gene][xdatatype]:
                    val= datasets[xdataset][gene][xdatatype][xdatakey]
                    if val not in taboo:
                        vals.append(val)
    return vals

def checkTaboo(pair):
    '''Given a pair of x, y values, return True if both values are allowed, False if one or both values are taboo'''
    
    x,y= pair[0], pair[1]
    if (x in taboo) or (y in taboo):
        return False
    else:
        return True
    
def getScatterVals(geneset='all'):
    '''
    return values for genes. Use this one if there is only one DE-Seq file to access and we're getting the baseMeans for FP, tot, or TE
    Allow homology dicitonary. key organism should be the first dataset.
    Implemented such that if multiple homologs found it will duplicate the x value and report all the y values (homologs)
    '''
    datasets= args['datasets']
    comparex= args['comparex']
    comparey= args['comparey']
   
    if args['attD']!=None:
        attD= unpickle(args['attD'])
    
    if args['y2h']!=None:
        homologyD= unpickle(args['y2h'])
    else:
        homologyD= None
        
    if homologyD==None:
        genespace= set(datasets[0].keys()) #only go through genes that are in all the dictionaries
        for i in range(1, len(datasets)):
            genespace=genespace.intersection(set(datasets[i].keys()))
    else:
        genespace= set(datasets[0].keys())
        
    if geneset=='all':
        geneset= genespace #this should be fine as long as all the datasets have the same set of genes
    
    if args['summaryfile']!=None:
        k= open(args['summaryfile'], 'w')
    else:
        k= open('test.txt', 'w')
        
    v1s=[]
    v2s=[]
    
    #figure out which values to collect for the x and yvalues:
    if comparex==['att']:
        comparex='att'
    else:
        xdataset= int(comparex[0])-1
        xdatatype= comparex[1]
        xdatakey= comparex[2]
    
    if comparey==['att']:
        comparey='att'
    else:
        ydataset= int(comparey[0])-1
        ydatatype= comparey[1]
        ydatakey= comparey[2]
    
    #collect the values            
    if (comparex=='att') or (comparey=='att'): #we're making the assumption that only one dataset can be an attribute here
        if comparex=='att':
            attkey=args['attKey']
        if comparey=='att':
            attkey=args['attKey']
        print 'attkey', attkey
        
        for gene in geneset:
            if gene in genespace: #Need to check b/c there's no guarantee that genes from the input geneset will be included in the data
                try:
                    v1= math.log(attD[gene][attkey],10)
                    v2= datasets[ydataset][gene][ydatatype][ydatakey]
                    if checkTaboo((v1, v2)):
                        line= '%s\t%s\t%s\t%s\n' % (gene, '', v1, v2)
                        k.write(line)
                        v1s.append(v1)
                        v2s.append(v2)
                except:
                    continue
        if comparey=='att': #reverse the x and y designations
            v1s= v2s.copy()
            v2s= v1s.copy()
        
    elif homologyD!=None:
        print 'finding homologs'
        k.write('%s\t%s\t%s\t%s\n' % ('yeastID', 'humanID', 'yeastval', 'humanVal'))
        for gene in geneset:
            
            if gene in genespace: #Need to check b/c there's no guarantee that genes from the input geneset will be included in the data
              
                try:
                    v1= datasets[xdataset][gene][xdatatype][xdatakey]
                    v2vals=[]
                    humanIDs= list(homologyD[gene])
                    for p in humanIDs:
                        v2vals.append(datasets[ydataset][p][ydatatype][ydatakey]) #this is a set, possibly containing >1 homolog
                    if len(v2vals)==1:
                        v1temp= [v1]
                        v2temp= [v2vals[0]]
                    else:
                        v1temp=[v1]*len(v2vals) #duplicate the x-value as many times as necessary
                        v2temp= v2vals #v2vals should be in the same order as their corresponding human ids
                        
                    for j in range(0, len(v1temp)):
                        if checkTaboo((v1temp[j], v2temp[j])):
                            line= '%s\t%s\t%s\t%s\n' % (gene, humanIDs[j], v1temp[j], v2temp[j])
                            k.write(line)
                            v1s.append(v1temp[j])
                            v2s.append(v2temp[j])
                except:
                    continue
        k.close()
    
    else:
        k.write('%s\t%s\t%s\t%s\n' % ('ID1', 'ID2', 'val1', 'val2'))
        
        for gene in geneset:
            if gene in genespace: #Need to check b/c there's no guarantee that genes from the input geneset will be included in the data
                
                try:
                    v1= datasets[xdataset][gene][xdatatype][xdatakey]
                    v2= datasets[ydataset][gene][ydatatype][ydatakey]
                    
                    if checkTaboo((v1, v2)):
                        line= '%s\t%s\t%s\t%s\n' % (gene, '', v1, v2)
                        k.write(line)
                        v1s.append(v1)
                        v2s.append(v2)
                except:
                    continue
    
    return v1s, v2s

def getMultiScatterVals(geneset='all'):
    '''
    return values for genes. Use this one if there is only one DE-Seq file to access and we're getting the baseMeans for FP, tot, or TE
    B/c of the way this is processed later, I only want to include genes with legitimate values in all of the libraries
    include in an array like this [[x1, x2...], [x1, ...], [], [],..]
    '''
    
    comparex= args['comparex']
    comparey= args['comparey']
    datasets= args['datasets']
   
    num=len(datasets)
    allVals=[[] for i in range(num)] #[[],[],[],[],..]
    
    genespace= set(datasets[0].keys()) #only go through genes that are in all the dictionaries
    for i in range(1, len(datasets)):
        genespace=genespace.intersection(set(datasets[i].keys()))
    if geneset=='all':
        geneset= genespace #this should be fine as long as all the datasets have the same set of genes
    
    #figure out which values to collect for the x and yvalues:
    xdatatype= comparex[1]
    xdatakey= comparex[2]
   
    for gene in geneset:
        if gene in genespace: #Need to check b/c there's no guarantee that genes from the input geneset will be included in the data
            tempvals=[] #temporary space to store all the vals
            for i in range(0, num):
                try:
                    v= datasets[i][gene][xdatatype][xdatakey]
                    
                    if v not in taboo:
                        tempvals.append(v)
                except:
                    continue
                
            if len(tempvals)==num: #All the datasets passed
                for i in range(0, num):
                    allVals[i].append(tempvals[i])
        
    return allVals

def transformIt(x='', y='', log=2):
    '''Log transform either one or two vectors'''
    
    if y=='': #only one dataset provided
        n1s=[]
        for i in range(0, len(x)):
            try:
                n1= math.log(x[i], log)
                n1s.append(n1)
            except:
                continue
        return n1s
        
    else:
            
        n1s=[]
        n2s=[]
    
        for i in range(0, len(x)):
            try:
                n1= math.log(x[i], log)
                n2= math.log(y[i], log)
            
                n1s.append(n1)
                n2s.append(n2)
            except:
                continue
        
        return n1s, n2s
    
def histifyCDF(vals):
    
    #num_bins= math.sqrt(len(vals))
    num_bins=100
    n_counts, bin_edges = np.histogram(np.array(vals), bins=num_bins, density=False)
    
    #make cdf
    cdf = np.cumsum(n_counts)
    scale = 1.0/cdf[-1] #scale 1/(total counts)
    ncdf = scale * cdf #normalized cdf , multiply everything by the scale factor, makes the max=1

    #make histogram:
    nhist= scale*n_counts
    
    return ncdf, nhist, bin_edges

def getMW(bg_vals, subset_vals):
    '''Determine the Mann Whitney U pvalue (one-sided test) for a change in ranks of bg and subset'''
    
    sig= stats.mannwhitneyu(bg_vals, subset_vals)
    pval=sig[1]
    #pval= math.log(sig[1],10)
    
    #since this operates on rank, need to figure out which list has higher ranks
    mut_rank=''
    combolist= copy.copy(bg_vals)
    combolist.extend(subset_vals)
    
    ranks=stats.rankdata(combolist)
    
    #figure out if the bg part of the list is < or > than the subset part of the list
    bgsum_avg= sum(ranks[0:len(bg_vals)])/len(bg_vals)
    subsetsum_avg= sum(ranks[len(bg_vals):])/len(subset_vals)
    
    #the average rank per gene in each list should indicate whether the bg or subset list is higher
    
    if bgsum_avg>subsetsum_avg: #this means that the bg values with have higher rank, i.e. generally be +ve compared to the subset values
        mut_rank='lower'
        down_pval= pval
        up_pval=1-pval
        return down_pval
        
    elif bgsum_avg<subsetsum_avg:
        mut_rank='higher'
        up_pval= pval
        down_pval=1-pval
        return up_pval
    else:
        mut_rank='NA'
        down_pval=''
        up_pval=''
        print 'equal ranks!'
        return None

def CDF(scatterobj=''):
    
    '''
    here dataType= 'FP', 'tot', 'TE'
    geneSetFiles= a list of geneset files, full path length required
    nesting for the data dictionaries is like this: gene->FP->'val'
    funcargs is a dictionary containing the arguments
    modify this so that we can give it a scatter object to use for histification if already collected elsewhere
    '''
    datasets= args['datasets']
    if args['filenum']== None:
        args['filenum']= int(args['comparex'][0])-1 #if we are not calling the list option, get directly
        
    if scatterobj=='':
        infinity= float('inf')
        taboo=set([infinity, -infinity, 'ND'])
        
        #define the options:
        geneSetFiles= args['genesets']
        legends=['all genes']
        if args['legends']!=None:
            legends.extend(args['legends'])
        mutants= args['mutants']
            
        #load gene names for each set:
        geneSets= addGenes(geneSetFiles)
            
        if args['transform']!=None:
            transform= args['transform']
            toTransform=True
            if transform=='log10':
                logbase=10
            elif transform=='log2':
                logbase=2
        else:
            transform=False
            toTransform=False
               
        #get gene counts for each group
        genecounts=[]
       
        #get bg vals first
        bg_vals=getVals()
        if toTransform==True:
            b1= transformIt(x=bg_vals, log=logbase)
            bg_vals=b1
        genecounts.append(len(bg_vals))
                          
        #get subset vals
        allsubsets=[]
        for subset in geneSets:
            subsetvals= getVals(geneset=subset)
            if toTransform==True:
                s1= transformIt(x=subsetvals, log=logbase)
                subsetvals=s1
            allsubsets.append(subsetvals)
            genecounts.append(len(subsetvals))
    
    else: #in this case we're turning a scatter plot into a cdf, so don't calculate
        genecounts= scatterobj.genecounts
        legends=['all genes']
        legends.extend(scatterobj.legends)
        mutants= scatterobj.mutants
        bg_vals= scatterobj.allVals[0]
        allsubsets= scatterobj.allVals[1:]
        
    #get histogram for each
    allCdfs=[]
    allHists=[]
    allEdges=[]
    allVals= []
    pvals=[]
    
    allVals.append(bg_vals)
    
    ncdf, nhist, bg_edges= histifyCDF(bg_vals)
    allCdfs.append(ncdf)
    allEdges.append(bg_edges)
    allHists.append(nhist)
    
    for subsetvals in allsubsets:
        if args['mw_pval']==True:
            pval=getMW(bg_vals, subsetvals)     
        else:
            pval=stats.ks_2samp(bg_vals,subsetvals)[1] #this is a two-sided p-value
        ncdf, nhist, s_edges= histifyCDF(subsetvals)  #since bin_edges are different for each dataset, have to save both
        allCdfs.append(ncdf)
        allHists.append(nhist)
        allEdges.append(s_edges)
        pvals.append(pval)
        allVals.append(subsetvals)
    
    #try making this one object and then pickle it:
    cdf= histObject(allCdfs, allHists, allEdges, pvals, legends, genecounts, allVals, mutants=mutants)
    cdf.otype='Histo'
        
    return cdf

def ratio(v1, v2):
    ''' return v2/v1 or 'ND' if v1=0'''
    try:
        r= float(v2)/v1
    except ZeroDivisionError:
        r= 'ND'
    return r
    
def findNewVals3(modify='comparex'):
    '''Given a function such as (ratio 1 TE A_r1 1 TE A_r2), perform the function (i.e. mean, ratio) on each gene
    and add to the datasets. Allow data from different datasets, which would be of the form (ratio 1 TE baseMeanA 2 TE baseMeanB)'''
    
    if args[modify][0].startswith('('):
        functions={'ratio':ratio}
        datasets= args['datasets']
        fields= args[modify]
        #clean up the fields by removing the outer parentheses:
        s= ' '.join(fields)
        s= s.strip('(').strip(')')
        fields= s.split(' ')
        
        functionname= fields[0]
        function= functions[functionname]
        
        datanum1= int(fields[1])-1 #the 0-based index of the dataset in datasets
        datatype1= fields[2]
        key1= fields[3]
        
        datanum2= int(fields[4])-1
        datatype2= fields[5]
        key2= fields[6]
        
        for gene in datasets[datanum1]:
            try:
                v1= datasets[datanum1][gene][datatype1][key1]
                v2= datasets[datanum2][gene][datatype2][key2]
                r= ratio(v1, v2)
                datasets[datanum1][gene][datatype1][functionname]= r
                
            except:
                continue
        
            args['datasets']= datasets #update datasets
            args[modify]=[fields[1], datatype1, functionname] #update comparex or comparey #how does it know which one to access??
    
    else:
        return
    
def Scatter():
    
    '''
    here dataType= 'FP', 'tot', 'TE'
    geneSetFiles= a list of geneset files, full path length required
    nesting for the data dictionaries is like this: gene->FP->'val'
    funcargs is a dictionary containing the arguments
    
    given two datasets, a datatype (FP, tot, TE), and one or more genesets (optional), return two vectors containing values
    all genes for which there is data
    '''
    findNewVals=False
   
    #define the options:
    geneSetFiles= args['genesets']
    legends= args['legends']
    comparex= args['comparex'] #determine if it's in parentheses--> this means some other function needs to be performed on it
    comparey= args['comparey']

    for b in ['comparex', 'comparey']:
        findNewVals3(modify=b)
   
    if args['y2h']!=None:
        homoloD= unpickle(args['y2h'])
    else:
        homoloD= None
        
    if args['summaryfile']!=None:
        outfile= args['summaryfile']
    else:
        outfile= None
        
    if args['transform']!=None:
        transform= args['transform']
        toTransform=True
        if transform=='log10':
            logbase=10
        elif transform=='log2':
            logbase=2
    else:
        transform=False
        toTransform=False
    
    allPlotData=[] #store the datasets as list of lists [[x1, y1], [x2, y2], ...]

    #load gene names for each set:
    geneSets= addGenes(geneSetFiles)    
    genecounts=[]
    
    v1, v2= getScatterVals()

    if toTransform==True:
        v1, v2= transformIt(x=v1, y=v2, log=logbase)
    allPlotData.append([v1, v2])
    genecounts.append(len(v1))
    
    for subset in geneSets:
        s1, s2= getScatterVals(geneset=subset)
        if toTransform==True:
            s1, s2= transformIt(x=s1, y=s2, log=logbase)
        allPlotData.append([s1, s2])
        genecounts.append(len(s1))
    #get pearson and spearman stats for each scatter:
    spear_rvalues=[]
    spear_pvalues=[]
    
    pear_rvalues=[]
    pear_pvalues=[]
    
    for comp in allPlotData:
        spear= stats.spearmanr(comp[0], comp[1])
        pear= stats.pearsonr(comp[0], comp[1])
        
        spear_rvalues.append(spear[0])
        spear_pvalues.append(spear[1])
        
        pear_rvalues.append(pear[0])
        pear_pvalues.append(pear[1])
        
        #first element is r value
        #second element is pvalue
    #make the scatter object
    scatter= scatterObject(allPlotData, spear_rvalues, spear_pvalues, pear_rvalues, pear_pvalues, legends=legends, genecounts=genecounts, mutants=args['mutants'], order=[])
    scatter.otype='Scattero'
    return scatter
    
    
def ScatterCDF():
    '''Use this to run the Scatter function to get the needed ratios, then run the getScatterRatio()
    function to convert it to a CDF object'''
    
    scatterobj= Scatter()
    cdfobj= getScatterRatio(scatterobj)
    cdfobj.otype= 'Histo'
    return cdfobj

def isGreaterThan(x, y):
    '''return True if y>x, else False'''
    if y>x:
        return True
    else:
        return False
    
def isLessThan(x, y):
    '''return True if y<x, else False'''
    if y<x:
        return True
    else:
        return False
    
def ScatterList():
    '''
    Similar to ScatterCDF, except return a list of the genes meeting the ScatterListObject specifications, rather than making an object
    Since I want the gene identities, can't simply use Scatter()
    '''
    
    operation={'greater':isGreaterThan, 'less':isLessThan}
    #define the options:
    legends= args['legends']
    comparex= args['comparex'] #determine if it's in parentheses--> this means some other function needs to be performed on it
    datanum= int(comparex[1])-1
    datatype= comparex[2]
    functionname= comparex[0].strip('(')
    co= args['scatterOutputList_co']
    direction= args['scatterOutputList_direction']
 
    #datasets[datanum1][gene][datatype1][functionname]= r
    passedGenes=set()
    findNewVals3(modify='comparex')
    for gene in args['datasets'][datanum]:
        if datatype in args['datasets'][datanum][gene]:
            if functionname in args['datasets'][datanum][gene][datatype]:
                v= args['datasets'][datanum][gene][datatype][functionname]
                if operation[direction](co, v):
                    passedGenes.add(gene)
    
    outname= '%s.txt' % args['outname']
    f=open(outname, 'w')
    for g in passedGenes:
        f.write('%s\n' % g)
    f.close()
    
def getOrder(num):
    '''return a list of paired indices to be plotted together'''
            
    plotPositions=[] #make this into a list of tuples of the form ([x,y], n) where x are the x values, y values and n is the spot in the grid
        
    range_y= range(0, num)
    range_x= range(0, num)
        
    counter=1
    for i in range_y:
        for j in range_x:
            pair= [i,j]
            if i==j:
                counter+= len(range_x)-j
                break
            else:
                plotPositions.append((pair, counter))
                counter+=1
                
    return plotPositions
    

def getScatterRatio(obj):
    '''If you create a scatter object but then want to get ratios between the different objects and then save them as a CDF object
    this will convert the scatter object into a CDF object'''
    
    xypairs= obj.allPlotData
    allVals=[]
    for i in range(0, len(xypairs)): #for each series
        xvals=[]
        pair=xypairs[i]
        for j in range(0, len(pair[0])):
            slope=pair[1][j]-pair[0][j] #log2(y)/log2(x)= log2(y-x)
            xvals.append(slope)
            
        allVals.append(xvals)
    
    obj.allVals= allVals
    cdfobj=CDF(obj) #convert this object into a CDF object
    return cdfobj


def writeScatterGenes(cdfobj):
    '''
    write the gene set that meets the cutoff requirement for the scatterOutputList
    '''
    
    vals= cdfobj.allVals[0] #only do this for the background vals for now
    
def MultiScatter():
    #define the options:
    geneSetFiles= args['genesets']
    legends= args['mutants']
    comparex= args['comparex']
    filenum= args['filenum']
    if args['transform']!=None:
        transform= args['transform']
        toTransform=True
        if transform=='log10':
            logbase=10
        elif transform=='log2':
            logbase=2
    else:
        transform=False
        toTransform=False
    
    #load gene names for each set:
    geneSets= addGenes(geneSetFiles)
    
    unpickledData= args['datasets'] #for the scatter() function, we will assume that there are only two datasets
    num= len(unpickledData)
    
    allVals= getMultiScatterVals() #removed support for geneset plotting in subsets
    order=getOrder(num)
    genecounts=[] #not implemented
    
    #GO THRU THE PAIRS AND GET THE R VALUES FOR EACH
    
    #get pearson and spearman stats for each scatter:
    spear_rvalues=[]
    spear_pvalues=[]
    
    pear_rvalues=[]
    pear_pvalues=[]
    
    
    for i in order:
        xindex= i[0][0]
        yindex= i[0][1]
        spear= stats.spearmanr(allVals[xindex], allVals[yindex])
        pear= stats.pearsonr(allVals[xindex], allVals[yindex])
                
        spear_rvalues.append(spear[0])
        spear_pvalues.append(spear[1])
        
        pear_rvalues.append(pear[0])
        pear_pvalues.append(pear[1])
        
    #make the scatter object
    scatter= scatterObject(allVals, spear_rvalues, spear_pvalues, pear_rvalues, pear_pvalues, legends, genecounts, order=order)
    scatter.otype= 'multiScattero'
    return scatter

def filterGenes():
    '''Given a minimum number of counts, remove the FP, tot, and TE as required if not enough counts in the libraries,
    My interpretation of Nick's binomial figure is that in order to get a reasonable estimate of the ratio, you need 128 counts total (between both conditions)
    therefore, for total, sum of all 4 total libraries (mutant and wt replicates)= 128
    '''
    co= args['counts']
    datasets= args['datasets']
    
    if args['filenum']!=None:
        tofilter=[args['filenum']]
    else:
        tofilter= range(0, len(datasets))
        
    datatypes=['FP', 'tot']
    if args['count_keys']!=None:
        datakeys= args['count_keys']
    else:
        datakeys=['A_counts_r1', 'A_counts_r2', 'B_counts_r1', 'B_counts_r2']
    #only allow length of datakeys to be 2 or 4
    #allow datakeys to be ['A_counts_r1', 'A_counts_r2' if we are plotting reproducibility]
    
    for i in tofilter:
        d= datasets[i]
        for gene in d:
            temp={}
            for datatype in datatypes:
                try: #datatype is already absent from the dictionary if there were no values to begin with
                    temp[datatype]=[]
                    total=0
                    for key in datakeys:
                        val= d[gene][datatype][key]
                        temp[datatype].append(val)
                        total+= val
                    
                    if total<co: #if total reads don't meet the cutoff, then delete this key
                        del d[gene][datatype]
                except:
                    continue
    
            #delete the TE entries unless both FP and total still exist
            if ('FP' in d[gene]) and ('tot' in d[gene]):
                #still have to check that A_r1, A_r2 for fp and total>= co    
                #test if both WT and mutant fp vs. total also pass the filter test:
                if len(datakeys)==4:
                    if ((sum(temp['FP'][0:2]) + sum(temp['tot'][0:2])) < co) or ((sum(temp['FP'][2:]) + sum(temp['tot'][2:])) < co):
                        del d[gene]['TE']
                        
                elif len(datakeys)==2:
                    if ((sum(temp['FP'][0:2]) + sum(temp['tot'][0:2])) < co):
                        del d[gene]['TE']
                    
            else:
                try:
                    del d[gene]['TE']
                except:
                    continue
        datasets[i]=d
        
    return datasets

def main(argList):
    
    global args
    
    objectList={'CDF':CDF, 'Scatter':Scatter, 'MultiScatter':MultiScatter, 'ScatterCDF':ScatterCDF, 'ScatterList':ScatterList}
    
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('object', help='type of object that you want to make, see objectList above')
    parser.add_argument('outname', help='outname of the file')
    parser.add_argument('-comparex', nargs='+', help='put the 1-based position of the library to get the vals, then datatype and the datakey, i.e. 1 FP val')
    parser.add_argument('-comparey', nargs='+', help='put the 1-based position of the library to get the vals, then datatype and the datakey, i.e. 1 FP val')
    parser.add_argument('-y2h', help='this should be a pickled dictionary mapping yeastID->set(homologous ensembl ids)')
    parser.add_argument('-summaryfile', help='name for a summary file')
    parser.add_argument('-genesets', nargs='+', help='list of genesets (as tab-delimited files) to be examined')
    parser.add_argument('-legends', nargs='+', help='labels for the genesets')
    parser.add_argument('-mutants', nargs='+', help='labels for the mutants, if multiple files used to make the object or object list')
    parser.add_argument('-list', action='store_true', help='make a list of objects, not a single object, looping in order through the input dictionaries')
    parser.add_argument('-transform', help='transform the data to logscale, support log2 and log10')
    parser.add_argument('-attD', help='enter an attribute dict mapping: gene->key->value, for example to look at transcript length vs. expression')
    parser.add_argument('-attKey', help='the name of the key for the values you want to collect in the attD')
    parser.add_argument('-filenum', type=int, help='this is the 0-based index of the filenumber to manipulated, generally called from a loop if automatedly doing a lot of mutants with the -list option')
    parser.add_argument('-infiles', nargs='+', help='the list of input files which are the pickled experiment dictionairies')
    parser.add_argument('-counts', type= int, default=0, help='this is the minimum number of counts that a gene must have in order to be included in the analysis')
    parser.add_argument('-count_keys', nargs='+', help='these are the keys the use for filtering, will override the defaults')
    parser.add_argument('-scatterOutputList_co', type=float, help='the value cutoff to make it in to the list')
    parser.add_argument('-scatterOutputList_direction', help='select one of greater or less to indicate if you want values to be greater or less than the cutoff')
    parser.add_argument('--mw_pval', action='store_true', help='calculate the Mann Whitney U pvalue for difference between groups, only implemented for CDF class')
    ar= parser.parse_args(args=argList)
    args= vars(ar)
    
    #add the unpickled datafiles to the args:
    args['datasets']=[]
    for i in args['infiles']:
        d= unpickle(i)
        args['datasets'].append(d)
        
    #remove genes from dictionary that do not meet the minimum counts cutoff:
    if args['counts']>0: #we want to filter the libraries for read counts    
        args['datasets']= filterGenes()
    
    obj= objectList[args['object']]()
    return obj, args['outname']
  
  
if __name__ == '__main__':
    main(sys.argv[1:])
