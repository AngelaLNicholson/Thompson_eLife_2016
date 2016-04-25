#150511 Mary K. Thompson
#Purpose: to make a plot of GO category ORF lengths vs. deltaTE to show that all short ones have decreased TE and we're not cherry-picking the mito groups
#GO dictionary format: ['P', 'C', 'F']->'0005291': {'genes': set(['YGR191W']), 'parents': set(['0005287', '0005290']), 'children': set([]), 'name': 'high affinity L-histidine transmembrane transporter activity'}

import sys
import cPickle as pickle

black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)
colors = [black, blue, reddishPurple, vermillion, orange, skyBlue, yellow, blue, yellow]
box_colors=[blue, reddishPurple, skyBlue, orange]
colorConvert={'black':black, 'orange':orange, 'skyBlue':skyBlue, 'bluishGreen':bluishGreen, 'yellow':yellow, 'blue':blue, 'vermillion':vermillion, 'reddishPurple':reddishPurple}

import numpy as np
import plotExtra
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import scipy.stats as stats

def arrangebyLen(treefile, lenfile, datafile, lenco):
    
    f= open(treefile, 'r')
    tree=pickle.load(f)
    f.close()
    
    g= open(lenfile, 'r')
    lens= pickle.load(g)
    g.close()
    
    h=open(datafile, 'r')
    d= pickle.load(h)
    h.close()
    
    bigD={}
    
    thisTree= tree['C']
    for gocat in thisTree:
        goname= thisTree[gocat]['name']
        theseLens=[]
        theseTEs=[]
        for gene in thisTree[gocat]['genes']:
            try:
                x1= lens[gene]['ORFLen']
                x2= d[gene]['TE']['val']
                theseLens.append(x1)
                theseTEs.append(x2)
            except:
                continue
        
        if len(theseLens)<lenco: #so if we set lenco=21, then 21 gene groups should pass
            continue
        else:
            bigD[goname]={}
            bigD[goname]['medianLen']= np.median(theseLens)
            bigD[goname]['medianTE']= 100*2**np.median(theseTEs)-100
        
    print 'num cats', len(bigD)
    return bigD

def plotD(d, outname):
    '''
    Plot the median deltaTEs of different GO categories arranged by median ORF length
    Can we make these so that now every GO category with mitochondria is red and ribosome is green
    Just plot over the current scatter plot, but need to keep track of the index so that we can do that
    '''
    
    l=[(d[name]['medianLen'], d[name]['medianTE'], name) for name in d]
    l.sort() #arrange by shortest to longest ORF length
    
    vals= [i[1] for i in l]
    names= ['%s (%s nt)' % (i[2], i[0]) for i in l]
    lens= [i[0] for i in l]
    
    #find the indices where the median lengths go > 500, 1000, 1500, 2000
    toFind=[500, 1000, 1500, 2000, 2500, 3000, 3500]
    p=0
    boundaryIs=[]
    for i in range(0, len(lens)):
        thisbound=toFind[p]
        if lens[i]>thisbound:
            boundaryIs.append(i+1) #add 1 to make 1-based
            p+=1
    
    newLabs= toFind[0:len(boundaryIs)]
    
    mitonames=[]
    mitovals=[]
    mitox=[]
    
    ribonames=[]
    ribovals=[]
    ribox=[]
    
    for j in range(0, len(names)):
        if 'mitochondria' in names[j]:
            mitonames.append(names[j])
            mitovals.append(vals[j])
            mitox.append(j+1) #since x is 1-based, then index is +1
        elif 'ribosom' in names[j]:
            ribonames.append(names[j])
            ribovals.append(vals[j])
            ribox.append(j+1) #since x is 1-based, then index is +1
        
    fig= plt.figure(figsize=(2.6,2))
    ax= fig.add_subplot(1,1,1)
    x= range(1, len(vals)+1)
 
    d1= ax.scatter(x, vals, s=5, color="k", edgecolors='none')
    d2= ax.scatter(mitox, mitovals, color=colorConvert["reddishPurple"], s=5, edgecolors='none')
    d3= ax.scatter(ribox, ribovals, color=colorConvert["blue"], s=5, edgecolors='none')
    
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_visible(False)
    
    ax.tick_params(axis='y', length=2, labelsize=5, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
    ax.tick_params(axis='x', length=2, labelsize=5, pad=2)
    ax.tick_params(top='off', bottom='on', right='off', left='on') #
    
    #plotExtra.addticks(ax,boundaryIs,newLabs,pos='x')
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel('median ORF length (nt)', fontsize=7)
    
    ax.set_xticks(boundaryIs)
    ax.set_xticklabels(newLabs)
    plt.legend([d1,d2,d3], ['all', 'mitochondrial', 'ribosomal'], loc=4, fontsize= 5, scatterpoints=1, markerscale=1.5)
    plt.axhline(y=0, linestyle='--', color='k', linewidth=0.5)
    plt.ylim(ymin=-25, ymax=10)
    
    plt.margins(0.02)
    plt.ylabel('percent change TE', fontsize=7)
    
    plt.savefig(outname + '.pdf', dpi=600, facecolor='w', edgecolor='w', format='pdf', transparent=True)
    
def main():
    
    treefile, lenfile, datafile, outname, lenco= sys.argv[1:]
    lenco= int(lenco) #gene group must have at least this many genes to be included
    d= arrangebyLen(treefile, lenfile, datafile, lenco)
    f=open('%s.p' % outname, 'w')
    pickle.dump(d, f)
    f.close()
    plotD(d, outname)

    
    '''
    treefile, lenfile, datafile, outname, lenco= sys.argv[1:]
    lenco= int(lenco) #gene group must have at least this many genes to be included
    d= arrangebyLen(treefile, lenfile, datafile, lenco)
    f=open('%s.p' % outname, 'w')
    pickle.dump(d, f)
    f.close()
    
    
    #for testing the plotting
    infile, outname= sys.argv[1:]
    f= open(infile, 'r')
    d= pickle.load(f)
    f.close()
    
    plotD(d, outname)
    '''
main()