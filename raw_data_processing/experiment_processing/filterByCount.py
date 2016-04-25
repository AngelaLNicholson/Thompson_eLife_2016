#140617 Mary K. Thompson
#purpose: filter the experiment objects by total read count so that these can easily be used for plotting, etc.

import sys
import cPickle as pickle
import os

def filterGenes(d, co):
    '''Given a minimum number of counts, remove the FP or total if not enough counts in each library,
    My interpretation of Nick Ingolia's binomial figure is that in order to get a reasonable estimate of the ratio, you need 128 counts total (between both conditions)
    therefore, for total, sum of all 4 total libraries (mutant and wt replicates)= 128
    '''
    datatypes=['FP', 'tot']
    datakeys=['A_counts_r1', 'A_counts_r2', 'B_counts_r1', 'B_counts_r2']
    
    for gene in d:
        temp={}
        for datatype in datatypes:
            if datatype not in d[gene]: #was already filtered out at an earlier step
                continue
            temp[datatype]=[]
            total=0
            for key in datakeys:
                val= d[gene][datatype][key]
                temp[datatype].append(val)
                total+= val
         
            if total<co: #if total reads don't meet the cutoff, then delete this key
                del d[gene][datatype]
        
        #delete the TE entries unless both FP and total still exist
        if ('FP' in d[gene]) and ('tot' in d[gene]):
            #still have to check that A_r1, A_r2 for fp and total>= co    
            #test if both WT and mutant fp vs. total also pass the filter test:
            if ((sum(temp['FP'][0:2]) + sum(temp['tot'][0:2])) < co) or ((sum(temp['FP'][2:]) + sum(temp['tot'][2:])) < co):
                del d[gene]['TE']
            
        else:
            try:
                del d[gene]['TE']
            except:
                continue
    return d

def report(d, co):
    '''Report how many genes passed cutoff in each category-- FP, total and TE'''
    
    datatypes=['FP', 'tot', 'TE']
    for datatype in datatypes:
        passedgenes=0
        for gene in d:
            if datatype in d[gene]:
                passedgenes+=1
        
        print '%s: %s genes passed with a cutoff of %s reads' % (datatype, passedgenes, co)
             

def main(argList):
    
    params= argList[0]
    prefix= os.path.join(params.outdir, params.bigoutname)
    
    for exp in range(0, len(params.reporder)):
	WT= params.WT[exp]
        for strain in params.reporder[exp]:
            if strain != WT:
                filename= '%s_%s_allData.p' % (prefix, strain)
                f= open(filename, 'r')
                d= pickle.load(f)
                f.close()
                
                d= filterGenes(d, params.co)
                report(d, params.co)
                
                newoutname= '%s_%s_allData_filtered.p' % (prefix, strain)
                g= open(newoutname, 'w')
                pickle.dump(d, g)
                g.close()

if __name__ == '__main__':
    main(sys.argv[1])