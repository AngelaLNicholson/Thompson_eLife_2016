#140627 MKT
# parse the .mp file output by bowtie1 after mapping to the rRNA index and perform downstream analyses
# need to input the strain name to do the analysis currently (sigma1, sigma2, human), since I already have the positions
# mapped for sigma, going to treat s288c the same way since I think the rRNA sequences should be very similar, but feel
# free to add an S288C version below

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cPickle as pickle
import os
import argparse
import cPickle as pickle
from collections import defaultdict
import numpy

#oligo sites based on the SGD-derived sequences for yeast
yeast_brar={'oligo2': {'start': 5753, 'length': 26, 'target': 'RDN37-1  '}, 'oligo3': {'start': 6419, 'length': 25, 'target': 'RDN37-1  '}, 'oligo1': {'start': 3989, 'length': 19, 'target': 'RDN37-1  '}}
yeast_eichorn={'Oligo 3': {'start': 6423, 'length': 21, 'target': 'RDN37-1  '}, 'Oligo 2': {'start': 5757, 'length': 21, 'target': 'RDN37-1  '}, 'Oligo 1': {'start': 3987, 'length': 21, 'target': 'RDN37-1  '}}
human_eichorn= {'Oligo 3': {'start': 4111, 'length': 21, 'target': 'rRNA28S'}, 'Oligo 2': {'start': 190, 'length': 21, 'target': 'rRNA28S'}, 'Oligo 1': {'start': 4748, 'length': 21, 'target': 'rRNA28S'}, 'Oligo 7': {'start': 4037, 'length': 21, 'target': 'rRNA28S'}, 'Oligo 6': {'start': 126, 'length': 21, 'target': 'rRNA5.8S'}, 'Oligo 5': {'start': 95, 'length': 21, 'target': 'rRNA5S'}, 'Oligo 4': {'start': 735, 'length': 21, 'target': 'rRNA18S'}}
celegans_doug={}

#store info about the rRNA targeting oligos
oligod={}
oligod['yeast']={}
oligod['yeast']['brar']= yeast_brar
oligod['yeast']['eichorn']= yeast_eichorn
oligod['human']={}
oligod['human']['eichorn']= human_eichorn
oligod['celegans']={}
oligod['celegans']['doug']= celegans_doug

#store info about the rRNA regions for each species
rRNAs={}
rRNAs['yeast']={'RDN37-1':{'length': 6845}, 'RDN5-1':{'length': 122}}
rRNAs['human']={'rRNA18S':{'length':1869}, 'rRNA28S':{'length':5070}, 'rRNA5S':{'length': 119}, 'rRNA5.8S':{'length':153}}
rRNAs['celegans']= {'ce6_rmsk_LSU-rRNA_Celrange_chrI_15060286-15061137':{'length':852}, 'ce6_rmsk_SSU-rRNA_Celrange_chrI_15062129-15063775':{'length':1647}, 'ce6_rmsk_LSU-rRNA_Celrange_chrI_15064825-15068333':{'length':3509},\
    'ce6_rmsk_SSU-rRNA_Celrange_chrI_15069326-15070972':{'length':1647}, 'ce6_rmsk_LSU-rRNA_Celrange_chrV_11938588-11938660':{'length':73}}

def plot_histogram(frag_dist, outName):
    
    f= open('%s_fragdist.txt' % outName, 'w')
    header='%s\t%s\t%s\n' % ('length', 'num_fragments', 'fraction')
    f.write(header)
    
    plt.rcParams['pdf.fonttype'] = 42

    x_labels=frag_dist.keys()
    x_labels.sort()
    bar_locations=numpy.arange(len(x_labels))
    
    width=.5
  
    vals=[]
    Total=0
    for bin in x_labels:
        Total+=frag_dist[bin]
    for bin in x_labels:
        percentage=frag_dist[bin]/float(Total)
        vals.append(percentage)
        f.write('%s\t%s\t%s\n' % (bin, frag_dist[bin], percentage))
  
    plot=plt.bar(bar_locations[:20], vals[:20], width, color='blue')
    plt.ylabel('Fraction of Reads')
    plt.xlabel('Read Length')
    plt.title(os.path.basename(outName).split('.')[0])
    plt.xticks(bar_locations[:20], x_labels[:20])
    plt.savefig('%s_fragdist.pdf' % outName, transparency=True, format='pdf')
    f.close()

def parseIn(inFile, outdir):
    
    frag_dist={} #keep track of number of reads of each length
    
    basename=os.path.basename(inFile).split('.')[0]
    outname= os.path.join(outdir, basename)
    RDNMappings={}
    RDNMappings['all']={}
    
    thisID=''    
    f= open(inFile, 'r')
    for line in f:
        
        fields= line.strip('\n').split('\t')
    
        ID= fields[0] #need to make this different for first one
        thisRegion=fields[2] #this will tell you 18S, 25S, etc.
        index= int(fields[3])
        seq= fields[4]
        read_length= len(seq)
        
        if read_length not in RDNMappings:
            RDNMappings[read_length]={}
            frag_dist[read_length]=0
            
        if thisRegion not in RDNMappings[read_length]:
            RDNMappings[read_length][thisRegion]=defaultdict(int)
            
        if thisRegion not in RDNMappings['all']:
            RDNMappings['all'][thisRegion]=defaultdict(int)
        
        RDNMappings[read_length][thisRegion][index]+=1
        RDNMappings['all'][thisRegion][index]+=1
            
        if ID != thisID:
            thisID=ID #we hit a new ID, count the reads
            frag_dist[read_length]+=1
            
    f.close()
    
    f=open('%s_allLens.p' % outname, 'w')
    pickle.dump(RDNMappings, f)
    f.close()
    return frag_dist, RDNMappings

def plotContam(Mappings, outName, oligod, rRNAd):
    '''
    Make a plot for each rRNA locus for each length group
    '''
    #find total number of all rRNA- mapping reads:
    totalrRNA=0
    for length in Mappings:
	for t in Mappings[length]:
	    totalrRNA+= sum(Mappings[length][t].values())
	
    for length in Mappings:
        #print 'length', length
        Alldensities={}
        
        for transcript in rRNAd:
            if transcript not in Mappings[length]: #it's possible that there will be no reads mapping to the transcript
                continue
            fileName= '%s_%s_%s' % (outName, length, transcript) 
            tlen= rRNAd[transcript]['length']
            Alldensities[transcript]={}
            Alldensities[transcript]['x']=range(0, tlen)
            Alldensities[transcript]['y']=[]
            Alldensities[transcript]['normalized']=[]
            Alldensities[transcript]['density']=[0.0]*tlen
    
            #get all the reads for each region
            for i in range(0, tlen):
                if i in Mappings[length][transcript]:
                    Alldensities[transcript]['y'].append(Mappings[length][transcript][i])
                else:
                    Alldensities[transcript]['y'].append(0)
            
            for i in range(0, tlen):
                Alldensities[transcript]['normalized'].append(float(Alldensities[transcript]['y'][i])*100/totalrRNA)
            
            #add estimated overlap from 3' bases, estimate with 25mer length reads
            if length=='all':
                addlen=28
            else:
                addlen= length
            for i in range(0, tlen):
                for j in range(addlen):
                    try:
                        Alldensities[transcript]['density'][i+j]+=Alldensities[transcript]['normalized'][i]
                    except:
                        break #we hit the end of the annotated region
            #plot the contaminants
            fig= plt.figure()
            
            ax= fig.add_subplot(1,1,1)     
            ax.plot(Alldensities[transcript]['x'], Alldensities[transcript]['density'], color='b')
            ax.fill_between(Alldensities[transcript]['x'], 0, Alldensities[transcript]['density'], facecolor='b')
            plt.ylim(ymin=0)
            plt.ylabel('percent of rRNA-mapping reads')
            plt.xlim(xmax= tlen-1)
            plt.title('%s contamination' % transcript)
           
            #now add the location of all the targeting oligos:
            for i in oligod:
                start= oligod[i]['start']
                end= start+oligod[i]['length']-1
               
                ax.axvline(x=start, color='r', linewidth=0.5, linestyle='--')
                ax.axvline(x=end, color='r', linewidth=0.5, linestyle='--')
        
            plt.savefig(fileName+'.pdf', dpi=600, facecolor='w', edgecolor='w', format='pdf')
            #plt.clf()        
            
def main(argList):
    
    parser= argparse.ArgumentParser(description= 'parse command line args')
    
    parser.add_argument('mpFile', help='the .mp file produced from mapping the rRNA reads')
    parser.add_argument('outputDir', help= 'the output directory where you want to put the results of this run')
    parser.add_argument('-species', help='enter human or yeast')
    parser.add_argument('-oligos', help='enter brar or eichorn')
    	
    ar=parser.parse_args(args=argList)
    args= vars(ar)
    
    if not os.path.exists(args['outputDir']):
        os.makedirs(args['outputDir'])
    
    outname= os.path.join(args['outputDir'], os.path.basename(args['mpFile']).strip('.mp'))
    
    #parse the .mp file, write an output pickle    
    frag_dist, RDNMappings= parseIn(args['mpFile'], args['outputDir'])
    
    #plot length of all the rRNA fragments
    plot_histogram(frag_dist, outname)
    
    #plot the rRNA contaminants of each length and overlay the rRNA targeting oligos from the selected species        
    theseOligos = oligod[args['species']][args['oligos']]
    thisrRNA= rRNAs[args['species']]
    
    plotContam(RDNMappings, outname, theseOligos, thisrRNA)
    
if __name__=='__main__':
    main(sys.argv[1:])