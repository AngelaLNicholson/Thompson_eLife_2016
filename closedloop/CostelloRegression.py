#150522 Mary K. Thompson
#purpose: Plot TE changes in mutants of mRNA groups from Costello et al., 2015. Also plot ORF lengths of those groups.

import sys
import cPickle as pickle
import csv
import numpy as np
import argparse

#import plottingFxns, assume that they are kept in the same directory structure as in github repo
moduleDir='../common/'
sys.path.append(moduleDir)

import linRegress as lr
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats
import matplotlib.ticker as plticker

RPfile="/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/Paper_planning/RACK1_paper/Figures/FIGUREDATA/F1/geneLists/RP_genes_PV.txt"
MRPfile="/Users/MK/Desktop/Gilbertlab/digital_lab_notebook/Paper_planning/RACK1_paper/Figures/FIGUREDATA/F1/geneLists/mitoribo_keyword.txt"
groupOrder=['Group 1', 'Group 2', 'Group 3A', 'Group 3B', 'Group 4A', 'Group 4B', 'Group 4C']

black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)

simple_colors = [reddishPurple, vermillion, black]
all_colors=[orange, bluishGreen, reddishPurple, yellow, vermillion, skyBlue, blue]

def cleanAxes(ax):
    '''Remove top and right lines from the graph'''
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    
    ax.tick_params(axis='y', length=2, labelsize=6, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
    ax.tick_params(axis='x', length=2, labelsize=6, pad=2)
    ax.tick_params(top='off', bottom='on', right='off', left='on') #
    
    return ax

def KaiserSmooth(x, beta):
    """ kaiser window smoothing
    from http://glowingpython.blogspot.co.uk/2012/02/convolution-with-numpy.html"""
    window_len=11
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(),s,mode='valid')
    
    return y[5:len(y)-5]

def parseList(infile):
    thisset=set()
    
    f= open(infile, 'r')
    for line in f:
        thisset.add(line.strip())
    f.close()
    return thisset

def readCostello(csvfile):
    '''
    Parse the Costello S5 table in csv format. Groups I and II have only one subdivision, III and IV have multiple
    '''
    
    f= open(csvfile, 'rU')
    reader= csv.reader(f)
    header=reader.next() #remove header
    
    GroupSets={0:{}, 1:{}, 2:{}, 3:{}}
    activeGroups=['Group 1', 'Group 2', 'Group 3', 'Group 4'] #name of groups that are being actively added to the dict now
    groupOrder=['Group 1', 'Group 2', 'Group 3A', 'Group 3B', 'Group 4A', 'Group 4B', 'Group 4C']
    
    rownum=0
    for row in reader:
        
        g1=row[0]
        g2=row[4]
        g3=row[8]
        g4=row[12]
        
        theseGenes=[g1, g2, g3, g4]
     
        for i in range(0, len(theseGenes)):
            
            if theseGenes[i]=='':
                if rownum==0: #on the first row if there is '' then this one does not have subgroups
                    GroupSets[i][activeGroups[i]]=set() #initialize these with Group 1, Group 2
            
            elif theseGenes[i].strip().upper().startswith('GROUP'):
                number= theseGenes[i].strip().split(' ')[1].upper() #reformat this to be the appropriate name
                for k in groupOrder:
                    if number in k:
                        cleanName=k
                        
                GroupSets[i][cleanName]= set()
                activeGroups[i]=cleanName
            
            else:    
                GroupSets[i][activeGroups[i]].add(theseGenes[i])
       
        rownum+=1

    f.close()
    
    #it makes sense to add by column, 0, 1, 2, 3 initially, but now I want to unnest them for easier manipulation later
    newGroups={}
    for k in GroupSets:
        for l in GroupSets[k]:
            newGroups[l]= GroupSets[k][l]
                
    return newGroups


def getGeneLists(order, cogroups):
    '''Get the sets of genes in each group'''
    neworder=[]
    for i in range(0, len(order)):
        geneset=set()
        for j in order[i]: #order[i] is a dict {name:set([group1, group2])}, j is the only key, order[i][j] should be the set([group1, group2])
            for p in order[i][j]:
                geneset= geneset.union(cogroups[p])

        #we want new order to mimic the order list, but with geneset instead of a list
        neworder.append({j:geneset}) #this should work b/c j is the only key
    return neworder

def convertToPercentage(allDs):
    '''convert all the log2 values to percentage from the nested list of lists structure'''
    for i in allDs:
        for j in allDs[i]['sens']: #this will be 'all', 'noRPs', 'noRPsorMRPs'
            newvs=[]
            for l in allDs[i]['sens'][j]:
                newl=[]
                for p in l:
                    v=100*2**p-100
                    newl.append(v)
                newvs.append(newl)
            allDs[i]['sens'][j]= newvs
                
    return allDs

def getAtts(cogroups, lenD, RPs, MRPs, d, outname, convertToPercent=False):
    '''
    Get the ORF lengths and asc1 sensitivity values for each group
    '''
    
    allDs={}
    keys=['all', 'simple']
    for k in keys:
        allDs[k]= {}
        for i in ['lens', 'sens']:
            allDs[k][i]=  {'all':[], 'noRPs':[], 'noRPsorMRPs':[]}
       
    groupOrder=[{'Group 1':set(['Group 1'])}, {'Group 2':set(['Group 2'])}, {'Group 3A':set(['Group 3A'])}, {'Group 3B':set(['Group 3B'])}, {'Group 4A':set(['Group 4A'])}, {'Group 4B':set(['Group 4B'])}, {'Group 4C':set(['Group 4C'])}]
    simpleOrder=[{'closed-loop enriched, repressor de-enriched':set(['Group 3A', 'Group 3B'])}, {'closed-loop enriched':set(['Group 4A'])}, {'other':set(['Group 1', 'Group 2', 'Group 4B', 'Group 4C'])}]
    
    #legends will simply be the associated list
    groupLegends=[i.keys()[0] for i in groupOrder]
    simpleLegends=[i.keys()[0] for i in simpleOrder]
    
    #get genesets for these groups:
    groupOrderGenes= getGeneLists(groupOrder, cogroups)
    simpleOrderGenes= getGeneLists(simpleOrder, cogroups)
    ordD={'all':groupOrderGenes, 'simple':simpleOrderGenes}
    
    #now want to go through and make a list of the values in that order for each category, all, noRPs, noRPsorMRPs
    for i in ['all', 'simple']:
        for j in ['all', 'noRPs', 'noRPsorMRPs']:
            thisOrder= ordD[i]
            for p in range(0, len(thisOrder)):
                
                ORFLens=[]
                asc1Sens=[]
                geneset= thisOrder[p].values()[0]
                
                for gene in geneset:
                    
                    if j=='all':
                        try:
                            x=float(lenD[gene]['ORFLen'])
                            y=float(d[gene]['TE']['val'])
                            ORFLens.append(x)
                            asc1Sens.append(y)
                        except KeyError:
                            continue
                        
                    if j=='noRPs':
                        if (gene in RPs):
                            continue
                        try:
                            x=float(lenD[gene]['ORFLen'])
                            y=float(d[gene]['TE']['val'])
                            ORFLens.append(x)
                            asc1Sens.append(y)
                        except KeyError:
                            continue
                        
                    if j=='noRPsorMRPs':
                        if (gene in RPs) or (gene in MRPs):
                            continue
                        try:
                            x=float(lenD[gene]['ORFLen'])
                            y=float(d[gene]['TE']['val'])
                            ORFLens.append(x)
                            asc1Sens.append(y)
                        except KeyError:
                            continue
                
                #append a list for each genelist
                allDs[i]['lens'][j].append(ORFLens) #this is where you put the orflens
                allDs[i]['sens'][j].append(asc1Sens)
    
    if convertToPercent==True:
        allDs=convertToPercentage(allDs)
    return allDs

def regressOut(allDs, how='total'):
    '''Regress ORF length out of the asc1-M1X sensitivity and see how adjusted asc1-M1X TE distribution looks.
    If how==byGroup, then each Group is regressed separately to get the residuals. Otherwise if how==total then
    regression done overall with all groups'''
        
    #regress groups separately
    if how=='byGroup':
        residuals=[]
        for i in range(0, len(allDs['simple']['lens']['all'])):
            a, b, RR, res= lr.linreg(allDs['simple']['lens']['all'][i], allDs['simple']['sens']['all'][i]) #res should give you the values of y with x regressed out
            residuals.append(res)
    else:
        #regress groups together, this is what is shown in the paper
        bigsensL=[]
        biglenL=[]
        lensoflists=[]
        for i in range(0, len(allDs['simple']['lens']['all'])):
            lensoflists.append(len(allDs['simple']['lens']['all'][i]))
            bigsensL.extend(allDs['simple']['sens']['all'][i])
            biglenL.extend(allDs['simple']['lens']['all'][i])
        
        a, b, RR, res= lr.linreg(biglenL, bigsensL) #res should give you the values of y with x regressed out
        
        #need to resplit into correct sizes
        residuals=[]
        starti=0
        for i in range(0, len(lensoflists)):
            endi=starti+lensoflists[i]
            residuals.append(res[starti:endi])
            starti= endi
                
    allDs['simple']['resis']= {}
    allDs['simple']['resis']['all']= residuals
        
    return allDs

def plotHist(hist, legends, outname, extraHists=(), smooth=True, linehistfill=True, beta=20, bin_range=(-2, 2, 0.05), x_range=(), plottype='sens'):
    '''
    hist is a list of lists containing values, legends is the name associated with them, including (n=x) number of genes
    extraHists is an optional tuple with hist lists and then the legend names assigned to them and then the indices of the hists to plot that will be plotted after
    the first set of hists in dotted lines of the same color
    '''
    fig = plt.figure(figsize=(2,2))
    ax1= fig.add_subplot(1,1,1)
    dataseries=[]
    
    if len(legends)>3:
        colors=all_colors
    else:
        colors=simple_colors
    
    allVals=[] #all values not split into lists
    for l in hist:
        for k in l:
            allVals.append(k)
  
    #now add the extraHists to the list:
    original_len= len(hist)
    if extraHists!=():
        for i in extraHists[2]: #add these indices only
            hist.append(extraHists[0][i])
        legends.extend(extraHists[1])
    
    #update legends with n= information
    new_legends=[]
    for i in range(0, len(legends)):
        new_legends.append('%s n=%s' % (legends[i], len(hist[i])))
    
    bins= np.arange(bin_range[0], bin_range[1], bin_range[2])
    
    allys=[]
    for i in range(0, len(hist)):
        ny, bin_edges = np.histogram(np.array(hist[i]), bins=bins, density=False)
        #use the middle of each xbin for plotting
        xs= [np.mean(bin_edges[m:m+2]) for m in range(0, len(bin_edges)-1)]
        
        totaly= float(sum(ny))
        ny_frac=[]
        for p in ny:
            frac= p/totaly
            ny_frac.append(frac)
        
        fillalphalevel=0.5
        linealphalevel=1
        
        #KAISER WINDOW SMOOTHING
        if smooth==True:
            yy=KaiserSmooth(ny_frac,beta)
        else:
            yy= ny_frac
        
        if i < original_len: #we're still in the original hist lengths
            s,=ax1.plot(xs, yy, alpha=linealphalevel, color=colors[i], linestyle='-', linewidth=0.5)
                 
            if linehistfill==True:
                    ax1.fill_between(xs, yy, alpha=fillalphalevel, color=colors[i])
                    testrect= plt.Rectangle((bin_edges[0],0), 0,0, facecolor=colors[i],alpha=fillalphalevel)
                    dataseries.append(testrect)
            else:
                    dataseries.append(s)
        else:
            s,=ax1.plot(xs, yy, alpha=linealphalevel, color=colors[extraHists[2][i-original_len]], linestyle='--', linewidth=0.5)
            dataseries.append(s)
        allys.append(yy)
        
    #find x limits based on where the data is:
    if x_range==():
        
        mindices=[]
        maxdices=[]
        for j in allys:
            for k in range(0, len(j)):
                
                if j[k]>0.01:
                    index=k
                    mindices.append(index)
                    break
        for j in allys:
            for k in range(len(j)-1, -1, -1):
                if j[k]>0.01:
                    index=k
                    maxdices.append(index)
                    break
        thismin=xs[min(mindices)]
        thismax=xs[max(maxdices)]
        ax1.set_xlim(xmin=thismin, xmax=thismax)
    else:
        ax1.set_xlim(xmin=x_range[0], xmax=x_range[1])

    ax1=cleanAxes(ax1)
    
    if plottype=='sens':
        loc1= plticker.MultipleLocator(base=20)
        ax1.axvline(x=0, linestyle='-', color='k', zorder=1, linewidth=0.5, alpha=0.3)
        plt.xlabel('deltaTE (asc1/WT)', fontsize=6)
    else:
        loc1= plticker.MultipleLocator(base=500)
        plt.xlabel('ORF length (nt)', fontsize=6)

    ax1.xaxis.set_major_locator(loc1)
    leg=ax1.legend(dataseries, new_legends, fontsize=6)
    plt.ylabel('fraction of genes', fontsize=6)
    plt.savefig(outname + '.pdf', dpi=600, facecolor='w', edgecolor='w', format='pdf', transparent=True)

def main(argList):
    parser= argparse.ArgumentParser()
    parser.add_argument('costello_file', help='CostelloS5.csv')
    parser.add_argument('lenpickle', help='pickled ORF and transcript lengths')
    parser.add_argument('datapickle', help='pickled experiment object, like either M1X_allData.p or 4G_allData.p')
    parser.add_argument('outname', help='outname for the graph')
    parser.add_argument('-how', default='total', help='method of regression, either total or byGroup')
    parser.add_argument('-beta', default=20, type=int, help='beta value for Kaiser window smoothing, lower=smoother')
    parser.add_argument('--plotlens', action='store_true', help='add if you want to plot the length distribution of the groups')
    parser.add_argument('--percentage', action='store_true', help='add if you want to convert from log2 fold changes to percentage')
    parser.add_argument('--mw', action='store_true', help='add if you want to report mann-whitney pvalues for differences between groups')
    
    ar= parser.parse_args(args=argList)
    args= vars(ar)
    
    cogroups=readCostello(args['costello_file'])
    f= open(args['lenpickle'], 'r')
    lenD= pickle.load(f)
    f.close()
    
    RPs= parseList(RPfile)
    MRPs= parseList(MRPfile)
    
    f= open(args['datapickle'], 'r')
    d= pickle.load(f)
    f.close()
    
    if args['percentage']==True:
        allDs= getAtts(cogroups, lenD, RPs, MRPs, d, args['outname'], convertToPercent=True)
        thisBinrange=(-50, 50, 2)
        thisXrange=(-40, 40)
    else:
        allDs= getAtts(cogroups, lenD, RPs, MRPs, d, args['outname'])
        thisBinrange=(-2, 2, 0.05)
        thisXrange=(-1, 1)

    allDs=regressOut(allDs, how=args['how'])
    
    #get pvalues for test of each group to the other group (only test simple closed loop and strong closed loop to the 'other' category)
    
    if args['mw']==True:
        for i in range(0, 2):
            normsig= stats.mannwhitneyu(allDs['simple']['sens']['all'][i], allDs['simple']['sens']['all'][-1])
            resisig= stats.mannwhitneyu(allDs['simple']['resis']['all'][i], allDs['simple']['resis']['all'][-1])
            orflensig= stats.mannwhitneyu(allDs['simple']['lens']['all'][i], allDs['simple']['lens']['all'][-1])
            print 'i', i
            print 'normsig', normsig
            print 'resisig', resisig
            print 'orflensig', orflensig
    '''
    #uncomment to get test values with no RPs or MRPs included
    if args['mw']==True:
        for i in range(0, 2):
            normsig= stats.mannwhitneyu(allDs['simple']['sens']['noRPsorMRPs'][i], allDs['simple']['sens']['noRPsorMRPs'][-1])
            resisig= stats.mannwhitneyu(allDs['simple']['resis']['noRPsorMRPs'][i], allDs['simple']['resis']['noRPsorMRPs'][-1])
            orflensig= stats.mannwhitneyu(allDs['simple']['lens']['noRPsorMRPs'][i], allDs['simple']['lens']['noRPsorMRPs'][-1])
            print 'i', i
            print 'normsig', normsig
            print 'resisig', resisig
            print 'orflensig', orflensig
    '''
    #full plot, graphs distribution for all groups from Costello paper as shown in supplemental data
    plotHist(allDs['all']['sens']['all'], groupOrder, ('%s_full' % args['outname']), smooth=True, bin_range=thisBinrange, x_range=thisXrange, beta=args['beta'], plottype='sens')

    #simple plot with residuals added
    plotHist(allDs['simple']['sens']['all'], ['cl+', 'cl', 'other'], ('%s_simple' % args['outname']), beta= args['beta'], extraHists= (allDs['simple']['resis']['all'], ['residual', 'residual','residual'], [0, 1, 2]), smooth=True, bin_range=thisBinrange, x_range=thisXrange, plottype='sens')
    
    if args['plotlens']==True:
        #plot all groups:
        plotHist(allDs['all']['lens']['all'], groupOrder, ('%s_alllens' % args['outname']), smooth=True, bin_range=(0, 5000, 100), x_range=(50,3700), plottype='lens')
        
        #plot simplified groups:
        plotHist(allDs['simple']['lens']['all'], ['cl+', 'cl', 'other'], ('%s_simplelens' % args['outname']), extraHists= (allDs['simple']['lens']['noRPsorMRPs'], ['no RPs', 'no RPs','no RPs'], [0, 1, 2]), smooth=True, bin_range=(0, 5000, 100), x_range=(50,3700), plottype='lens')    
    
if __name__ == '__main__':    
    main(sys.argv[1:])