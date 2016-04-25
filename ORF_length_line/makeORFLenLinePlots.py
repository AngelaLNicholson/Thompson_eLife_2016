#150827 Mary K. Thompson
#purpose: plot deltaTE, FP or total vs. ORF length for a given mutant

import sys
import cPickle as pickle
import numpy as np
import scipy.stats as stats
import csv
from collections import defaultdict
from scipy.stats import gaussian_kde
from scipy.stats import scoreatpercentile
from numpy import arange
import matplotlib.ticker as plticker
import argparse
import math
import copy

import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

RPgenefile='/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/Paper_planning/RACK1_paper/Figures/FIGUREDATA/F1/geneLists/RP_genes_PV.txt'
MRPgenefile='/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/Paper_planning/RACK1_paper/Figures/FIGUREDATA/F1/geneLists/mitoribo_keyword.txt'

black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)

colorConvert={'black':(0,0,0), 'orange':(230/255.0,159/255.0,0),'skyBlue':(86/255.0,180/255.0,233/255.0), 'bluishGreen':(0,158/255.0,115/255.0),\
    'yellow':(240/255.0,228/255.0,66/255.0), 'blue':(0,114/255.0,178/255.0), 'vermillion':(213/255.0,94/255.0,0), 'reddishPurple':(204/255.0,121/255.0,167/255.0)}

def parseGeneList(infile):
    
    s= set()
    f= open(infile, 'r')
    for line in f:
        gene=line.strip('\n')
        s.add(gene)
    return s

def splitByORFLenEqualGenesPlus(x, y, numgenes=100):
    '''Bin by ORF length according to number of genes contained in each, return averages of those intervals
    Now instead of saving the mean x value, I'm going to space them out on the x-axis by fraction of the genes quantified
    This is almost like dividing equally, but maybe slightly different on the ends because there might be a group on the end
    that isn't filled'''
    
    l=zip(x,y)
    l.sort()
    
    newx=[] #keep track of the new x values
    newy=[]
    newystd=[]
    
    bi=0 #keep track of bin index
    ngenes=0 #keep track of number of genes
    runningx=[] #running tally of gene lengths, average whenever we get to next bin and add to newx
    runningy=[]
    for i in range(0, len(l)):
        if ngenes/numgenes>bi:
            bi= ngenes/numgenes #new bin, add values for last bin, will currently only add last bin with a full 100 genes
            newx.append(np.mean(runningx))
            newy.append(np.mean(runningy))
            newystd.append(np.std(runningy))
            
            runningx=[]
            runningy=[]
            
        runningx.append([l[i][0]])
        runningy.append([l[i][1]])
        
        ngenes+=1
    
    num_genes=len(l)
    thisx= [((i+1)*100.0)/num_genes for i in range(0, len(newx))]

    return thisx, newy, newystd, newx #return the xvalue means as well so that we can mark on the x-axis the point at which the mean ORF length exceeds that value

def convertToPercentSimple(data, logbase):
    '''Simple convert to percentage points, return one list'''
    
    transformeddata=[]
    for d in data:
        v= (logbase**d)*100-100
        transformeddata.append(v)
    return transformeddata

def plotLine_general(datasets, legends, ORFLenD, colors, datakey, outname, exclude_ASC1=False, exclude_RPs=False):
    '''
    instead of plotting scatterplot, plot average of points at regular bins, here give
    the script a list of datasets of form d[gene]['TE']['val']=log2_deltaTE and an orf length dictionary,
    d[gene]['ORFLen']= orflength. legends= a list of legends to label the plots
    '''
    
    num_ysets= len(datasets)
    
    xvals=[]
    ysets=[]
    for i in range(num_ysets):
	ysets.append([])
	xvals.append([])
    
    #get all overlapping genes as the geneset
    genesets=[] #make list of genesets to iterate through when collecting data
    if exclude_RPs==False:
	allsets=[]
	for i in range(0, len(datasets)):
	    theseGenes=set()
	    for gene in datasets[i]:
		try:
		    if type(datasets[i][gene][datakey]['val'])==float:
			theseGenes.add(gene)
		except KeyError:
		    continue
		
	    allsets.append(theseGenes)
	
	newset=allsets[0]
	for i in range(1, len(allsets)):
	    newset=newset.intersection(allsets[i])
	
	
	if exclude_ASC1==True:
	    if 'YMR116C' in newset:    
		newset.remove('YMR116C')
	
	'''
	#add this temporarily to remove the eIF4G isoforms:
	if exclude_ASC1==True:
	    if 'YGL049C' in newset:
		newset.remove('YGL049C')
	    if 'YGR162W' in newset:
		newset.remove('YGR162W')
	'''	
	genesets.append(newset)
    
    else:
	for i in range(0, len(datasets)):    
	    newset=set(datasets[i].keys())
	    if exclude_ASC1==True:
		if 'YMR116C' in newset:
		    newset.remove('YMR116C')
	    genesets.append(newset)
	    
    #this adds them all separately now so that datasets might include slightly different genesets
    for i in range(0, len(ysets)):
	if exclude_RPs==True:
	    thisset= genesets[i]
	else:
	    thisset= genesets[0] #all the same
	    
	for gene in thisset:
	    try:
		x1= ORFLenD[gene]['ORFLen']
		y1= datasets[i][gene][datakey]['val']
	
		xvals[i].append(x1)
		ysets[i].append(y1)
	    except:
		continue
		
    spearman_vals=[]
    for j in range(0, len(ysets)):
	corr= stats.spearmanr(xvals[j], ysets[j])
	print 'corr', corr
	spearman_vals.append(corr)
    
    fig= plt.figure(figsize=(3,3))
    ax= fig.add_subplot(1,1,1)

    allx=[]
    ally=[]
    allstd=[]
    
    theseXmeans=[]
    for k in range(0, len(ysets)):
	x, yvals, ystd, xmeans= splitByORFLenEqualGenesPlus(xvals[k], ysets[k])
	theseXmeans.append(xmeans)
    
	allx.append(x)
	ally.append(yvals)
	allstd.append(ystd)
    
    xvals2plot=[]
    vals2plot=[]
    stds2plot=[]
    
    for i in range(0, len(ally)):
	yc= convertToPercentSimple(ally[i], 2)
	ystd= convertToPercentSimple(allstd[i], 2)
	xvals2plot.append(np.array(allx[i]))
	vals2plot.append(np.array(yc))
	stds2plot.append(np.array(ystd))
      
    series=[]
    for k in range(0, len(vals2plot)):
	ax.fill_between(xvals2plot[k], vals2plot[k]-stds2plot[k], vals2plot[k]+stds2plot[k], facecolor=colorConvert[colors[k]], edgecolor=colorConvert[colors[k]], alpha=0.15) #really doesn't seem to put this on the bottom??
	d=ax.scatter(xvals2plot[k], vals2plot[k], color=colorConvert[colors[k]], linewidth=0.5, s=5)

	series.append(d)
    
    ax.axhline(y=0, linestyle='-', color='k', zorder=1, linewidth=0.5, alpha=0.3)
    
    #still want to note ORF lengths with a second axis/set of ticks
    ax.set_ylim(-40, 40)
    
    #now that xmins can be slightly different b/c of different genesets, set as max of them
    firstxvals=[p[0] for p in xvals2plot]
    xmini= max(firstxvals)
    ax.set_xlim(xmin=xmini, xmax=1)

    ax.set_xlabel('fraction of genes', fontsize= 7)
    ax.set_ylabel('percent change %s' % datakey, fontsize=7)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
     
    newax= fig.add_axes(ax.get_position(), frameon=False)
    newax.tick_params(top='off', bottom='on', right='off', left='off')
    newax.set_yticks([])
    newax.set_xticks([])
    newax.set_xticklabels([])
    
    #set now with indices corresponding to the different ORF lengths
    toFind=[500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500]
    
    #find the indices that correspond to different ORF lengths
    xmeans=theseXmeans[0] #set to first one so that if we're plotting something accessory pts like noRPs, it will reference the original
    boundaryIs=[]
    newLabs=[]
    p=0
    for i in range(0, len(xmeans)): #note that for this new version of the plot, xmeans needs to be the same between different genesets
        thisbound=toFind[p]
        if xmeans[i]>thisbound:
	    boundaryIs.append(xvals2plot[0][i])
            p+=1
    
    newLabs= toFind[0:len(boundaryIs)]
    newax.set_xticks(boundaryIs)
    newax.set_xticklabels(newLabs)
    
    newax.get_xaxis().set_tick_params(which='both', direction='out')
    newax.tick_params(axis='x', length=2, labelsize=5, pad=2)
    
    loc1= plticker.MultipleLocator(base=20)
    loc2= plticker.MultipleLocator(base=0.2)

    ax.yaxis.set_major_locator(loc1)
    ax.xaxis.set_major_locator(loc2)

    ax.tick_params(axis='y', length=2, labelsize=5, pad=2) #format tick length and distance to labels
    ax.tick_params(axis='x', length=2, labelsize=5, pad=2)
    ax.tick_params(top='off', bottom='on', right='off', left='on') #

    #add the spearman correlation values to the legend
    for i in range(0, len(legends)):
	pval= spearman_vals[i][0]
	legends[i]= '%s r=%1.2f' % (legends[i], pval)
	    
    leg=ax.legend(series, legends, fontsize=5, scatterpoints=1)
    leg.draw_frame(False)
    
    plt.savefig('%s_line.pdf' % outname)
    
def excludeRPs(dataset, excludedSet):
    
    '''
    return the dataset with the genes in the excludedList removed
    '''
    
    #make new data dict excluding all the ones in the excluded set
    
    newd={}
    for i in dataset:
	if i in excludedSet:
	    continue
	else:
	    newd[i]= dataset[i]
    
    return newd
    
def main(argList):
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('outname')
    parser.add_argument('ORFLenD', help='this is the precompiled dictionary with ORFlens')
    parser.add_argument('-datafiles', nargs='+', help= 'add pickled data files, i.e. allData.p or allData_filtered.p')
    parser.add_argument('-legends', nargs='+', help='add strings in same order as datafiles to add as legends')
    parser.add_argument('-colors', nargs='+', help='name and order of the colors to plot')
    parser.add_argument('--exclude_RPs', action='store_true', help='add this if you want to exclude the RPs from the length correlation, will be applied to each datafile and plotted with unexcluded')
    parser.add_argument('--exclude_ASC1', action='store_true', help='add if you want to plot without ASC1 included in the plot')
    parser.add_argument('-datatypeKey', help='add FP, TE, total for which type of data you want to plot')	
    ar=parser.parse_args(args=argList)
    args= vars(ar)
    
    outname= args['outname']
    
    f= open(args['ORFLenD'], 'r')
    ORFLenD= pickle.load(f)
    f.close()
    
    datasets=[]
    for i in args['datafiles']:
	f= open(i, 'r')
	d= pickle.load(f)
	f.close()
	datasets.append(d)
    
    legends= args['legends'] #so for the case where we exclude RPs, len(datafiles)!= len(legends) at initial input

    excludedDs=[]
    if args['exclude_RPs']==True:
	RPs= parseGeneList(RPgenefile)
	MRPs= parseGeneList(MRPgenefile)
	allRPs= RPs.union(MRPs)
	
	excludedDs=[]
	#now generate a version of the dataset with those genes excluded:
	for i in datasets:
	    excludedD= excludeRPs(i, allRPs)
	    
	excludedDs.append(excludedD)
    
    
    if excludedDs!=[]:
	data2plot=[]
	for j in range(0, len(datasets)):
	    data2plot.append(datasets[j])
	    data2plot.append(excludedDs[j])
    else:
	data2plot= datasets
	
    plotLine_general(data2plot, legends, ORFLenD, args['colors'], args['datatypeKey'], outname, exclude_ASC1= args['exclude_ASC1'], exclude_RPs=args['exclude_RPs'])

    
if __name__ == '__main__':
    main(sys.argv[1:])
