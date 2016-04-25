#141028 Mary K. Thompson
#plot translational efficiency vs. mRNA attributes

import sys
import csv
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42

import scipy.stats as stats

def main():
    infile, outname= sys.argv[1:]
    
    atts= pd.DataFrame.from_csv(infile)
    
    theseCats= {'ORF length':'ORFLen', '5p UTR length':'median5', '3p UTR length':'median3', 'TE':'WT_TE', 'mRNA expression':'WT_rpkm',\
            '5p folding energy':'5pMFE', '3p folding energy': '3pMFE', 'polyA tail length': 'pALe', 'tRNA adapation Index':'tAI'}
    theseLimits={'mRNA expression':{'min':1, 'max':5},'TE':{'min':-1, 'max':1}, 'tRNA adapation Index':{'min':0, 'max':0.8}, '5p UTR length':{'min':0, 'max':400},\
                 '3p UTR length':{'min':0, 'max':400},'ORF length':{'min':0, 'max':6000}, '5p folding energy':{'min':-100, 'max':0}, '3p folding energy':{'min':-100, 'max':0},\
                 'polyA tail length':{'min':0, 'max':80}, 'ORF length2':{'min':0, 'max':6000}}
    toPlot=['mRNA expression','TE', 'tRNA adapation Index', '5p UTR length','3p UTR length','ORF length', '5p folding energy', '3p folding energy', 'polyA tail length']    
   
    j=1
    letters=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    fig= plt.figure(figsize=(8,8))
    for i in toPlot:
        ax= fig.add_subplot(3,3,j)
	
	mydata= atts[['deltaTEM1X', theseCats[i]]].dropna(how="any")
	vals= mydata.values
        
        d=ax.scatter(vals[:,0], vals[:,1], alpha=0.1)

	spearman, pvalue= stats.spearmanr(vals[:,0], vals[:,1])
	
        #set plot parameters
        ax.set_ylim(ymin=theseLimits[i]['min'], ymax=theseLimits[i]['max'])
        ax.set_xlim(xmin=-2, xmax=2)
        ax.set_ylabel('%s' % i, labelpad=2, fontsize=7)
        if j>6:
            ax.set_xlabel('%s' % 'log2_ TE (asc1-M1X/WT)', labelpad=2, fontsize=7)
        
        ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.spines['right'].set_visible(False)
        ax.text(0.8,0.85,'r=%1.3f\np=%1.3e' % (spearman, pvalue), transform= ax.transAxes, fontsize=5)

        ax.text(0.05,1, letters[j-1], transform= ax.transAxes, fontsize=9, horizontalalignment='left', verticalalignment='top')

        ax.tick_params(axis='y', length=0.8, labelsize=5, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
	ax.tick_params(axis='x', length=0.8, labelsize=5, pad=2)
        ax.tick_params(top='off', bottom='on', right='off', left='on') #
        j+=1
    
    fig.tight_layout()
    plt.savefig('%s.pdf' % outname, format='pdf', transparent=True)
    
    #plot part 2d, TE vs. ORF length
    fig= plt.figure(figsize=(2,2))
    ax= fig.add_subplot(1,1,1)
    
    mydata= atts[['ORFLen', 'WT_TE']].dropna(how="any")
    vals= mydata.values
    d= ax.scatter(vals[:,0], vals[:,1], s=0.2, color='k')
    
    #get the spearman correlations
    spearman, pvalue= stats.spearmanr(vals[:,0], vals[:,1])
  
    #set plot parameters
    ax.set_ylim(ymin=-1, ymax=1)
    ax.set_xlim(xmin=0, xmax=6000)
    ax.set_ylabel('translational efficiency (log10 TE)', labelpad=2, fontsize=7)
    ax.set_xlabel('ORF length (nt)', labelpad=2, fontsize=7)
    
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    ax.text(0.8,0.85,'r=%1.3f\np=%1.3e' % (spearman, pvalue), transform= ax.transAxes, fontsize=5)
    
    ax.tick_params(axis='y', length=2, labelsize=5, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
    ax.tick_params(axis='x', length=2, labelsize=5, pad=2)
    ax.tick_params(top='off', bottom='on', right='off', left='on') #
    
    plt.savefig('%s_TE_v_ORFLen.pdf' % outname, format='pdf', transparent=True)

main()