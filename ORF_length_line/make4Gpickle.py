#150523 Mary K. Thompson
#purpose: parse the eIF4G data from Park et al. 2011 and convert it to an pickled an experiment file with same structure as others

import sys
import math
import csv
import cPickle as pickle
import numpy as np

def fgMeanChange(M1, M2, M3, W1, W2, W3):
    
    try:
        val= np.mean([(float(M1)-float(W1)),(float(M2)-float(W2)), (float(M3)-float(W3))])
    except:
        val=''
        
    return val
    
def build4Gdata(fGfile):
    
    d={}
    f= open(fGfile, 'rU')
    reader= csv.reader(f)
    header= next(reader)
    for row in reader:
        genename= row[1]
        d[genename]={}
        for i in range(3, len(row)):
            d[genename][header[i]]= row[i]
           
    f.close()
    
    newd={}
    #calculate the deltaFP, deltaTot and deltaTE here:
    for gene in d:
        newd[gene]={}
        newd[gene]['FP']={}
        newd[gene]['tot']={}
        newd[gene]['TE']={}
        
        newd[gene]['FP']['val']= fgMeanChange(d[gene]['log4G-HP_I_14564002_532'], d[gene]['log4G-HP_II'], d[gene]['log4G-HP_III'], d[gene]['logWT-HP-I_14601302_532'],d[gene]['logWT-HP_II'],d[gene]['logWT-HP_III'])
        newd[gene]['tot']['val']= fgMeanChange(d[gene]['log4G-T-I_14592702_532'], d[gene]['log4G-T_II'], d[gene]['log4G-T_III'], d[gene]['logWT-T-I_14612202_532'], d[gene]['logWT-T_II'], d[gene]['logWT-T_III'])
        try:
            newd[gene]['TE']['val']= float(newd[gene]['FP']['val'])- float(newd[gene]['tot']['val'])

        except:
            newd[gene]['TE']['val']=''
  
    return newd

def build4GdataLP(fGfile):
    
    d={}
    f= open(fGfile, 'rU')
    reader= csv.reader(f)
    header= next(reader)
    for row in reader:
        genename= row[1]
        d[genename]={}
        for i in range(2, len(row)):
            d[genename][header[i]]= row[i]
           
    f.close()
    
    newd={}
    #calculate the deltaFP, deltaTot and deltaTE here:
    for gene in d:
        newd[gene]={}
        newd[gene]['TE']={}
	try:
            newd[gene]['TE']['val']= math.log(float(d[gene]['mean TE4g_LP/mean TEwt_LP']),2)
	except:
	    newd[gene]['TE']['val']=''
    return newd

def main():
    fG_HPfile, fG_LPfile, outname=sys.argv[1:]
    
    fgHP= build4Gdata(fG_HPfile)
    fgLP= build4GdataLP(fG_LPfile)

    hpfile= open('%s_hp.p' % outname, 'w')
    lpfile= open('%s_lp.p' % outname, 'w')
    
    pickle.dump(fgHP, hpfile)
    pickle.dump(fgLP, lpfile)
    
main()