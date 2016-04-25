#150822 Mary K. Thompson
#Purpose: To look at the relationship between ORF length, 5p, 3p UTR length, and deltaTE in asc1
#also going to plot ORF+3'UTR or ORF+5'UTR or ORF alone vs. deltaTE in M1X

import sys
import csv
import cPickle as pickle
import scipy.stats as stats
from operator import itemgetter

#import plottingFxns, assume that they are kept in the same directory structure as in github repo
moduleDir='../common/'
sys.path.append(moduleDir)

import pandas as pd
import analysesRunner
import graphPlotters
import numpy as np
import os
import linRegress as lr
import math

def unpickle(infile):
    f= open(infile, 'r')
    d= pickle.load(f)
    f.close()
    return d

def plotScatterPD(new_atts, ORFx, tcombox, y, outname, limit=None, set_bounds='wtTE'):
    '''
    given ORF lengths, UTR lengths, and deltaTE, correlate both with deltaTE and then perform partial correlations.
    The PD version, I'm going to pass it the names of the pandas columns
    '''
    
    mydata= new_atts[[ORFx, tcombox, y]]
    vals= mydata.values
    
    ORFx= vals[:,0]
    tcombox= vals[:,1]
    y= vals[:,2]
        
    if limit!=None: #if there is an ORF length limit, then collect the values < this limit:
        l= zip(ORFx, tcombox, y)
        l.sort()
        
        newl=[i for i in l if i[0]<=limit]
        
        ORFx= [i[0] for i in newl]
        tcombox= [i[1] for i in newl]
        y= [i[2] for i in newl]
    
    bounds={'ORF_length':{'y':[-2, 2], 'x':[0, 6000]}, 'combined_ORF_and_UTR':{'y':[-2, 2], 'x':[0, 6000]}, 'versus':{'y':[0, 6000], 'x':[0, 6000]},\
        'control_ORFLen':{'y':[-2,2], 'x':[-2000,2000]}, 'control_combo':{'y':[-2,2], 'x':[-2000,2000]}, 'control_TE':{'x':[0,6000], 'y':[0,6000]}}
    
    WTTEbounds={'ORF_length':{'y':[-1.5, 1.5], 'x':[0, 5000]}, 'combined_ORF_and_UTR':{'y':[-1.5, 1.5], 'x':[0, 5000]}, 'versus':{'y':[0, 5000], 'x':[0, 5000]},\
       'control_ORFLen':{'y':[-1.5,1.5], 'x':[-400,1000]}, 'control_combo':{'y':[-1,1], 'x':[-1000,400]}, 'control_TE':{'x':[0,5000], 'y':[0,5000]}}
    
    if set_bounds=='wtTE': #rescale y-axis to [-1, 1]
	bounds=WTTEbounds
        
    labels={'ORF_length':{'y':'log2 deltaTE', 'x':'ORF length (nt)'}, 'combined_ORF_and_UTR':{'y':'log2 deltaTE', 'x':'combined length (nt)'}, 'versus':{'y':'combined length (nt)', 'x':'ORF length (nt)'},\
        'control_ORFLen':{'y':'residual (deltaTE vs. ORF len)', 'x':'residual (combo len vs. ORF len)'}, 'control_combo':{'y':'residual (delta TE vs. combo len)', 'x': 'residual (ORF len vs. combo len)'}, 'control_TE':{'y':'residual (combo len vs. TE)', 'x':'residual (ORF len vs. TE)'}}
      
    toPlot={'ORF_length':[ORFx, y], 'combined_ORF_and_UTR':[tcombox, y], 'versus':[ORFx, tcombox]}
    
    thisOrder=['ORF_length', 'combined_ORF_and_UTR', 'versus']
    
    for i in thisOrder:
        X= toPlot[i][0]
        Y= toPlot[i][1]
	
        #Can modify args from another file like this:        
        graphPlotters.args['xmin']=bounds[i]['x'][0]
        graphPlotters.args['xmax']= bounds[i]['x'][1]
        graphPlotters.args['ymin']=bounds[i]['y'][0]
        graphPlotters.args['ymax']= bounds[i]['y'][1]
        graphPlotters.args['xlabel']= labels[i]['x']
        graphPlotters.args['ylabel']= labels[i]['y']
        graphPlotters.args['scattersize']=0.2
	graphPlotters.args['legendfont']=6
	graphPlotters.args['tickfont']=5
	graphPlotters.args['autolimits']=False
	graphPlotters.args['hidelegenddots']=True
	
        thisscatter= analysesRunner.scatterObject([[X,Y]], [stats.spearmanr(X,Y)[0]], [stats.spearmanr(X,Y)[1]], [stats.pearsonr(X,Y)[0]], [], genecounts=[len(X)])
        plot= graphPlotters.ScatterPlotter()
        plot.draw_subplot(thisscatter, (1,1,1))
	
	#override default pearson r value plotting, plot spearman r value
	plot.legends=['All (r=%1.2f)\np=%1.2e' % (thisscatter.spear_rvalues[0], thisscatter.spear_pvalues[0])]
        graphPlotters.formatIt(thisscatter, plot, 'Scatter')
        plot.save_plots('%s_%s_scatter' % (outname,i))
    
    #now get and plot partial correlations:
    toRegress={'control_ORFLen':[[ORFx, tcombox], [ORFx, y]], 'control_combo':[[tcombox, ORFx], [tcombox, y]], 'control_TE':[[y, ORFx],[y, tcombox]]}

    for i in toRegress:
        x1= toRegress[i][0][0]
        y1= toRegress[i][0][1]
        x2= toRegress[i][1][0]
        y2= toRegress[i][1][1]
        
        a, b, RR, res1= lr.linreg(x1, y1)
        a, b, RR, res2= lr.linreg(x2, y2)
        
        graphPlotters.args['xmin']=bounds[i]['x'][0]
        graphPlotters.args['xmax']= bounds[i]['x'][1]
        graphPlotters.args['ymin']=bounds[i]['y'][0]
        graphPlotters.args['ymax']= bounds[i]['y'][1]
        graphPlotters.args['xlabel']= labels[i]['x']
        graphPlotters.args['ylabel']= labels[i]['y']
        graphPlotters.args['scattersize']=1
        
        thisscatter= analysesRunner.scatterObject([[res1, res2]], [stats.spearmanr(res1,res2)[0]], [stats.spearmanr(res1,res2)[1]], [stats.pearsonr(res1,res2)[0]], [], genecounts=[len(res1)])
        plot= graphPlotters.ScatterPlotter()
        plot.draw_subplot(thisscatter, (1,1,1))
        graphPlotters.formatIt(thisscatter, plot, 'Scatter')
        plot.save_plots('%s_%s_scatter' % (outname,i))
        		        

def main():
    attcsv, lenpickle, outname= sys.argv[1:]
    lend= unpickle(lenpickle)
    
    if len(outname.split('/'))>1:
	thisPath=os.path.dirname(outname)
	if not os.path.exists(thisPath):
	    os.makedirs(thisPath)
	
    atts= pd.DataFrame.from_csv(attcsv)
    
    #1) Look at ORF length vs. txt length for the new data:
    #make this dataframe for getting ORF length or txt length vs. WT TE
    wt_atts= atts.copy()
    
    #get and add the txt lengths from here:
    txtlend={}
    for gene in lend:
	try:
	    txtlend[gene]= lend[gene]['median']
	except:
	    continue
    
    txt_len= pd.Series(txtlend)
    wt_result= pd.concat([wt_atts, txt_len], axis=1)
    wt_result= wt_result.dropna(subset=['ORFLen', 'WT_TE', 0])
    
    #Plot scatter plots:
    plotScatterPD(wt_result, 'ORFLen', 0, 'WT_TE', '%s_txtcomp_WT_TE' % outname)
    
    #2) Look at the effect of adding 3' UTR or 5' UTR to the ORF on deltaTE in M1X
    new_atts= atts.copy()
    new_atts['txt3']= new_atts['median3']+ new_atts['ORFLen'] #get ORF len + 3' side
    new_atts['txt5']= new_atts['median5']+ new_atts['ORFLen'] #get ORF len + 5' side
    
    new_result= pd.concat([new_atts, txt_len], axis=1)
    new_result= new_result.dropna(subset=['ORFLen', 0, 'deltaTEM1X', 'txt3', 'txt5'])
    
    #Plot scatter plots:
    plotScatterPD(new_result, 'ORFLen', 0, 'deltaTEM1X', '%s_txtcomp_deltaTE' % outname)
    plotScatterPD(new_result, 'ORFLen', 'txt3', 'deltaTEM1X', '%s_3comp_deltaTE' % outname)
    plotScatterPD(new_result, 'ORFLen', 'txt5', 'deltaTEM1X', '%s_5comp_deltaTE' % outname)

if __name__ == '__main__':
    main()
    
'''
#the lend looks like this:
{'average': 781.653164556962, 'median': 772.0, 'ORFLen': 597, 'sepLen': 777.0}

#Pelechano table S3 data looks like this:
chr	strand	gene	type	mTIF_number	median5	median3	sd5	sd3
1	-	YAL067C	Verified	1	33	222	NA	NA

'''