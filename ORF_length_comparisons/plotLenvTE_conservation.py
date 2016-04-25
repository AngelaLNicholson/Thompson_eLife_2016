#150606 Mary K. Thompson
#purpose: plot scatter plots of TE vs. ORF length in other organisms

import sys
import cPickle as pickle
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats

def plotScatter(X, Y, name, outname, species):
    
    bounds={'yeast':{'ORF length':{'y':[-1, 1], 'x':[0, 5000]}},'human':{'ORF length':{'y':[-1, 1], 'x':[0,6000]}}, 'mouse':{'ORF length':{'y':[-1, 1], 'x':[0,6000]}}, 'worm':{'ORF length':{'y':[-1, 1], 'x':[0,6000]}}}
    
    labels= {'yeast':{'ORF length':{'y':'log10 TE', 'x':'ORF length (nt)'}},'human':{'ORF length':{'y':'log10 TE', 'x':'ORF length (nt)'}}, 'mouse':{'ORF length':{'y':'log10 TE', 'x':'ORF length (nt)'}}, 'worm':{'ORF length':{'y':'log10 TE', 'x':'ORF length (nt)'}}}
    
    fig= plt.figure(figsize=(2,2))
    ax= fig.add_subplot(1,1,1)
    
    spearman, pvalue= stats.spearmanr(X, Y)
    
    d= ax.scatter(X, Y, s=0.2, color='k')
    
    thesebounds= bounds[species][name]
    theselabels= labels[species][name]
    
    ax.set_ylim(ymin=thesebounds['y'][0], ymax= thesebounds['y'][1])
    ax.set_xlim(xmin=thesebounds['x'][0], xmax= thesebounds['x'][1])
   
    ax.set_ylabel(theselabels['y'], labelpad=2, fontsize=7)
    ax.set_xlabel(theselabels['x'], labelpad=2, fontsize=7)
    
    ax.text(0.8,0.85,'r=%1.3f\np=%1.3e' % (spearman, pvalue), transform= ax.transAxes, fontsize=6)
    
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    
    ax.tick_params(axis='y', length=2, labelsize=5, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
    ax.tick_params(axis='x', length=2, labelsize=5, pad=2)
    ax.tick_params(top='off', bottom='on', right='off', left='on') #

    plt.savefig('%s.pdf' % (outname), format='pdf', transparent=True)
    plt.clf()


def runScatter(d, outname, species='yeast'):
    '''
    collect data, get correlations and pass to scatter plot function
    '''
    x=[]
    y=[]
    for gene in d:
       x.append(d[gene]['ORFLen'])
       y.append(d[gene]['TE'])
    
    plotScatter(x, y, 'ORF length', outname, species)
    
def main():
    pickledData, outname, thisSpecies= sys.argv[1:]
    
    f= open(pickledData, 'r')
    d= pickle.load(f)
    f.close()
        
    runScatter(d, outname, species=thisSpecies)
    
main()