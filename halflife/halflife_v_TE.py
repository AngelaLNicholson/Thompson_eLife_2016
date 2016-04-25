#150824 Mary K. Thompson
#purpose: make plots for TE vs. halflife

import sys
import csv
import cPickle as pickle
import scipy.stats as stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
plt.rcParams['pdf.fonttype'] = 42

import math
import numpy as np

def parseMiller(halffile):
    '''Parse Miller 2011 halflives'''
    d={}
    f= open(halffile, 'rU')
    f.readline()
    for line in f:
        row=line.strip('\n').split('\t')
        sysname= row[0]
        pAhalf= float(row[2])
        #totalhalf= float(row[3])
        d[sysname]={}
        d[sysname]['pA']=pAhalf
        #d[sysname]['tot']=totalhalf
        
    f.close()
    return d


def parseNey(halffile):
    '''Parse Neymotin 2014 halflives'''
    d={}
    f= open(halffile, 'rU')
    reader= csv.reader(f)
    header= reader.next()
    for row in reader:
        try:
            sysname= row[0]
            half= float(row[9])
            d[sysname]={}
            d[sysname]['pA']=half
        except:
            continue

    f.close()
    
    
    halflives= [d[gene]['pA'] for gene in d]
    print 'median', np.median(halflives)
    print 'average', np.mean(halflives)
    return d


def parsePres(halffile):
    '''Parse Presnyak 2015 halflives'''
    d={}
    f= open(halffile, 'rU')
    reader= csv.reader(f)
    header= reader.next()
    for row in reader:
        sysname= row[0]
        pAhalf= float(row[2])
        totalhalf= float(row[3])
        d[sysname]={}
        d[sysname]['pA']=pAhalf
        d[sysname]['tot']=totalhalf
        
    f.close()
    return d

def compareTEs(halfd, TEd, outname):
     
    x=[]
    y=[]
        
    for gene in halfd:
        try:
            x1= TEd[gene]
            y1= math.log(halfd[gene]['pA'], 10)
            x.append(x1)
            y.append(y1)
    
        except:
            continue
    
    corr=stats.spearmanr(x, y)
   
    print 'corr', corr
    
    fig= plt.figure(figsize=(2,2))
    ax= fig.add_subplot(1,1,1)
    d1=ax.scatter(x, y, s=0.2, color='k')
   
    ax.set_xlabel('TE')
    ax.set_ylabel('log10 half-life (minutes)')
    ax.set_ylim(ymin=0, ymax=3)
   
    ax.text(0.85, 0.85,'r=%1.3f' % corr[0], transform=ax.transAxes, fontsize=6)
            
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    
    ax.tick_params(axis='y', length=2, labelsize=5, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
    ax.tick_params(axis='x', length=2, labelsize=5, pad=2)
    ax.tick_params(top='off', bottom='on', right='off', left='on') #

    plt.savefig('%s_%s' % (outname,'halflife_TE.pdf'), format='pdf')

    plt.clf()
    
def readAttributes(attfile):
    
    d={}
    f= open(attfile, 'rU')
    reader= csv.reader(f)
    header= reader.next()
    for row in reader:
       
        try:
            WT_TE= float(row[6])
            thisID= row[0]
            d[thisID]=WT_TE
        except:
            continue
        
    return d

def main():
    halfMiller, halfPres, halfNey, attributefile, outname= sys.argv[1:]
    
    #Get halflife data:
    mill=parseMiller(halfMiller)
    pres=parsePres(halfPres)
    ney= parseNey(halfNey)
    
    #Get TE data:
    TEd= readAttributes(attributefile)
    
    compareTEs(mill, TEd, '%s_mill' % outname)
    compareTEs(pres, TEd, '%s_pres' % outname)
    compareTEs(ney, TEd, '%s_ney' % outname)

main()