#131007 Mary K. Thompson
#purpose: to run ks tests, hypergeometric, or mann-whitney test (as specified) on all of the GO category deltas, also write a file with the fold changes

import os
import cPickle as pickle
import scipy.stats as stats
import sys
import numpy
import math
import argparse
import operator
import copy

convertNamespace={'cellular_component':'C','molecular_function':'F','biological_process':'P'}


def benjaminiHochbergCorrection(pValDict):
    '''
    takes a dictionary mapping key to p value
    returns a dictionary of Benjamini-Hochberg corrected Q values

    Q = p * n / k, where n is the # of observations, and k is the rank of the particular p-value among all p-values
    '''
    qValues = {}
    sorted_p = sorted(pValDict.iteritems(), key=operator.itemgetter(1)) #this is going to sort low->high, therefore the low pvalue will have low rank
    sorted_p.reverse()
    
    n = len(sorted_p)
    for i in range(n):
        k = i+1
        q = sorted_p[i][1] * n / k
        qValues[sorted_p[i][0]] = q
    return qValues

def bonferroniCorrection(pValDict):
    '''
    takes a dictionary mapping key to p value
    returns a dictionary of Bonferroni corrected Q values

    Q = p * n, where n is the # of observations
    '''
    qValues = {}
    sorted_p = sorted(pValDict.iteritems(), key=operator.itemgetter(1))
    n = len(sorted_p)
    for i in range(n):
        q = sorted_p[i][1] * n
        qValues[sorted_p[i][0]] = q
    return qValues

def signPvals(d, ID):
    
    '''
    Make this now such that it only changes the sign of the pvalues to artificially make things with increases positive
    '''
    
    foldchange_down=d[ID]['down']['foldchange']
    foldchange_up= d[ID]['up']['foldchange']
    
    if foldchange<0:
        d[ID][s]['pval']= math.log(d[ID][s]['pval'], 10)
        d[ID][s]['benjamini']= math.log(d[ID][s]['benjamini'],10)
        d[ID][s]['bonferroni']= math.log(d[ID][s]['bonferroni'],10)
    elif foldchange>0:    
        d[ID][s]['pval']= -math.log(d[ID][s]['pval'], 10)
        d[ID][s]['benjamini']= -math.log(d[ID][s]['benjamini'],10)
        d[ID][s]['bonferroni']= -math.log(d[ID][s]['bonferroni'],10)
    
    d[ID][s]['foldchange']= foldchange
    
    return d, s

def writeOut(allDs, stanzas, category, outname):
    '''
    write an output file with all the pvalues. Make them positive if the median is increasing, negative if median is decreasing.
    '''
    
    header='\t'.join(['description', 'GO_ID', 'pval', 'benj', 'bonf', 'log2foldchange', 'num_enriched', 'num_total'])+'\n'

    for s in ['up', 'down']:
        outfilename= '%s_%s.txt' % (outname,s)
        g=open(outfilename, 'w')
        g.write(header)
            
        for ID in allDs:
            GOIDName='GO:%s' % ID
            description= stanzas[category][ID]['name']
            #for the MannWhitney test the 'up and down' should be perefectly symmetrical
            line=[description, GOIDName, '%1.3f' % math.log(allDs[ID][s]['pval'],10), '%1.3f' % math.log(allDs[ID][s]['benjamini'],10), '%1.3f' % math.log(allDs[ID][s]['bonferroni'],10), '%1.3f' % allDs[ID][s]['foldchange'], str(allDs[ID][s]['subset_enriched']),  str(allDs[ID][s]['subset_total'])]
    
            strline='\t'.join(line)+'\n'
            g.write(strline)
                
        g.close()
    
def countGenes(d, geneset, comparison, checkkey, coval):
    '''
    Check whether each gene meets the requirement
    '''
    up=0
    down=0
    total=0
    vals=[] #store all vals so can use to get the median, etc.
 
    for gene in geneset:
        if (gene in d) and (comparison in d[gene]):
            testval=abs(d[gene][comparison][checkkey])
            v= d[gene][comparison]['val']

            if testval >= coval: #if there is no change, then it won't be included
                if v>0:
                    up+=1
                elif v<0:
                    down+=1
                    
            vals.append(v)
            total+=1
     
    return up, down, total, vals

def testGOs(stanzas, dataname, d, outfile, pval_co=None, fold_change_co=None, min_genes=2, max_genes=10000, stat='ks'):
    '''
    Go through each comparison (FP, tot, TE) and each namespace, collect genes for each list and test enrichment
    change this now so the we can only assign genes to bg and subset by pval or foldchange, not both simultaneously
    '''
    
    if pval_co!=None:
        coval= pval_co
        checkkey= 'pval'
    else:
        coval= fold_change_co
        checkkey= 'val'
    
    j=open(stanzas, 'r')
    stanzas=pickle.load(j)
    j.close()
        
    for comparison in ['FP', 'tot', 'TE']:
        for category in stanzas: #like MF, CC, BP
            outname='%s_%s_%s_%s' %(outfile, dataname, comparison, category)
            
            allDs={}
            allGOs= stanzas[category].keys()
           
            allgenes= set(d.keys())
           
            bg_up, bg_down, bg_num, bg_vals= countGenes(d, allgenes, comparison, checkkey, coval)
            bg_median= numpy.median(bg_vals)
                  
            #Now test each GO category against the background set  
            for GO_ID in stanzas[category]:
                genes= list(stanzas[category][GO_ID]['genes'])
                num_genes= len(genes)
                subset_up, subset_down, subset_num, subset_vals= countGenes(d, genes, comparison, checkkey, coval)
                subset_median= numpy.median(bg_vals)
                
                if not min_genes<=subset_num<=max_genes: #check that it is between the gene number cutoffs
                    continue
                
                mean_fold_change= numpy.mean(subset_vals) #this is the GEOMETRIC MEAN
                subset_median= numpy.median(subset_vals) #why
                
                #Note here that we are not forcing different datasets to only include only the same set of genes, they will include all genes for which they have data--- perhaps we should change this
                if stat=='ks':
                    sig= stats.ks_2samp(bg_vals, subset_vals)
                    pval=sig[1]

                    if subset_median < bg_median:
                        median_is='lower'
                        down_pval= pval
                    elif subset_median > bg_median:
                        median_is='higher'
                        up_pval= pval
                    else:
                        median_is='ND'
                
                elif stat=='mw': #find the mann whitney u test pvalue
                    sig= stats.mannwhitneyu(bg_vals, subset_vals)
                    pval=sig[1]
                    
                    #since this operates on rank, need to figure out which list has higher ranks
                    mut_rank=''
                    combolist= copy.copy(bg_vals)
                    combolist.extend(subset_vals)
                    
                    ranks=stats.rankdata(combolist)
                    
                    #figure out if the bg part of the list is < or > than the subset part of the list
                    bgsum_avg= sum(ranks[0:len(bg_vals)])/len(bg_vals)
                    subsetsum_avg= sum(ranks[len(bg_vals):])/len(subset_vals)
                    
                    #the average rank per gene in each list should indicate whether the bg or subset list is higher
                    
                    if bgsum_avg>subsetsum_avg: #this means that the bg values with have higher rank, i.e. generally be +ve compared to the subset values
                        mut_rank='lower'
                        down_pval= pval
                        up_pval=1-pval
                        
                    elif bgsum_avg<subsetsum_avg:
                        mut_rank='higher'
                        up_pval= pval
                        down_pval=1-pval #is this OK? I think so, but not sure...
                        
                    else:
                        mut_rank='NA'
                        down_pval=''
                        up_pval=''
                        print 'equal ranks!'

                elif stat=='hg': #use Fischer's exact test to get the pvalue
                    
                    #note that here the bg does not include the subset, think this is the right way to do it.
                    a=subset_up
                    c=subset_num-subset_up
                    b=bg_up-subset_up #the number uniquely up in the bg
                    e=(bg_num-subset_num)-b #the number not up in the bg
                
                    up_mat=[[a,b],[c,e]] #used e instead of d b/c I already assigned d to the data dictionary
                    up_sig= stats.fisher_exact(up_mat, alternative='greater') #note, earlier versions of stats have bug in scipy.stats.fischer_exact
                   
                    
                    #to test downreg gene enrichment
                    a=subset_down
                    c=subset_num-subset_down
                    b=bg_down-subset_down
                    e=(bg_num-subset_num)-b
                 
                    down_mat=[[a,b],[c,e]]
                     
                    down_sig= stats.fisher_exact(down_mat, alternative='greater')
                    #since I only want to test for overrepresentation, use alternative=
                
                    up_pval=up_sig[1]
                    down_pval=down_sig[1]#log transform later
                    
                #Note here that I'm only choosing the best one (up or downregulated), although it is possible for this test to get a significant value for both (if bimodal distribution)
                allDs[GO_ID]={}
                allDs[GO_ID]['up']={}
                allDs[GO_ID]['up']['pval']= up_pval
                allDs[GO_ID]['up']['foldchange']=mean_fold_change
                allDs[GO_ID]['up']['subset_enriched']= subset_up
                allDs[GO_ID]['up']['subset_total']=subset_num
                
                allDs[GO_ID]['down']={}
                allDs[GO_ID]['down']['pval']= down_pval
                allDs[GO_ID]['down']['foldchange']=mean_fold_change
                allDs[GO_ID]['down']['subset_enriched']=subset_down
                allDs[GO_ID]['down']['subset_total']=subset_num
            
            #calculate the corrected pvalues:
            for s in ['up', 'down']:
                pValDict=dict(zip([k for k in allDs], [allDs[k][s]['pval'] for k in allDs]))
                benjamini_corrected= benjaminiHochbergCorrection(pValDict)
                bonfer_corrected=bonferroniCorrection(pValDict)
                for GO_ID in allDs:
                    allDs[GO_ID][s]['benjamini']= benjamini_corrected[GO_ID]
                    allDs[GO_ID][s]['bonferroni']= bonfer_corrected[GO_ID]
            
            if args['write_numgenes']==True:
                writeOut(allDs, stanzas, category, outname)
            
def main(argList):
    
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('outname', help= 'outname')
    parser.add_argument('pickledStanzas', help= 'a pickle of GO-annotated genes, as produced by parseObo.py')
    parser.add_argument('-files', nargs='+', help='list of pickled data files corresponding to mutants for analysis')
    parser.add_argument('-min_genes', type=int, default=3, help='minumum number of genes a GO category must have for inclusion, default=2')
    parser.add_argument('-max_genes', type=int, default=6000, help='maximum number of genes a GO category can have for inclusion, default= unlimited')
    parser.add_argument('-pval_co', type=float, help='a pvalue cutoff, in linear, not logscale to apply. If pval>pvalco for all samples, that GO category will not be written to the file')
    parser.add_argument('-fold_change_co', default=0, type= float, help= 'a fold change cutoff in logscale. If abs(foldchange)< fold_change_co for all samples, that GO category will not be written to the output file')
    parser.add_argument('--ks', action='store_true', help='get pvals from the 2 sample KS test, this is the defult behavior')
    parser.add_argument('--hg', action='store_true', help='get pvals from hypergeometric test')
    parser.add_argument('--mw', action='store_true', help='get pvals from the mann whitney u test')
    parser.add_argument('--write_numgenes', action='store_true', help='also write the number of genes enriched in each category and the total number of genes in that category. Only use if only doing one mutant')
    parser.add_argument('--make_gsea_in', action='store_true', help='write the gsea .gmx file for GSEA analysis')
   
    global args
    
    ar=parser.parse_args(args=argList)
    args= vars(ar)
        
    outname= args['outname']
    pickledStanzas= args['pickledStanzas']
    ds=[] #make this a list of data dictionairies
    datanames=[]
    infiles= args['files']
    for i in infiles:
        f= open(i, 'r')
        d= pickle.load(f)
        f.close()
        ds.append(d)
        name=os.path.basename(i).split('.')[0]
        datanames.append(name)

    for j in range(0, len(ds)): #runs them one-by-one, no overhead integration of the results
        datadict= ds[j]
        dataname= datanames[j]
        if args['ks']==True:
            testGOs(pickledStanzas, dataname, datadict, outname, pval_co=args['pval_co'], fold_change_co=args['fold_change_co'], min_genes=args['min_genes'], max_genes=args['max_genes'])
            
        if args['hg']==True:
            testGOs(pickledStanzas, dataname, datadict, outname, pval_co=args['pval_co'], fold_change_co=args['fold_change_co'], min_genes=args['min_genes'], max_genes=args['max_genes'], stat='hg')
        
        if args['mw']==True:
            testGOs(pickledStanzas, dataname, datadict, outname, pval_co=args['pval_co'], fold_change_co=args['fold_change_co'], min_genes=args['min_genes'], max_genes=args['max_genes'], stat='mw')
        

if __name__ == '__main__':
    main(sys.argv[1:])