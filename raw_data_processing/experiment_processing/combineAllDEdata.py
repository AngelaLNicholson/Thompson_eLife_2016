#140611 Mary K. Thompson
#purpose: to calculate translation efficiencies, changes between experimental samples, and combine data into experiment files (i.e. asc1_allData.p) to store the data

import sys
import cPickle as pickle
import csv
import os
import numpy as np
import math

taboo= set(['Inf', '-Inf'])

def stripLine(line, splitter):
    '''
    returns the components of a line with all whitespaces removed and split on splitter
    '''
    l=[]
    fields= line.split(splitter)
    for i in fields:
        l.append(i.strip())
    
    return l

def parseSfactors(sizefile):
    '''Get the sfactors from a DESeq run'''
    
    sfactors={}
    
    f= open(sizefile, 'r')
    header=f.readline() #nothing useful in the first line
    for line in f:
        libname, sf= line.strip('\n').split('\t')
        sfactors[libname]=float(sf)
    return sfactors


def parseDESeq(DEFile, countd, pc_normedd, repkeys, includeco=1):
    '''
    parse a DESeq output file and get the log2 values. If the log2 value=NA, there were no reads in any of the libraries.
    If the log2value=+/-inf, add the pseudocount values to all libraries with 0 reads, EXCEPT for cases where there is
    only 1 read that leads to the +/- inf value as this causes a non-sensical result
    '''
    
    DESeq_d={}
    
    f=open(DEFile, 'r')
    header=f.readline()
        
    for line in f:
	fields=line.strip('\n').split('\t')
	log2change= fields[6]
	if log2change!= 'NA':
	    if log2change not in taboo: #there was not a + or -ve infinity value, or zero in both libraries
		ID= fields[1]
		baseMeanA= float(fields[3])
		baseMeanB= float(fields[4])
		foldchange= float(fields[5])
		log2change= float(fields[6])
		pval= float(fields[7])
		padj= float(fields[8])
	
		DESeq_d[ID]={}
		DESeq_d[ID]['baseMeanA']= baseMeanA
		DESeq_d[ID]['baseMeanB']= baseMeanB
		DESeq_d[ID]['foldchange']= foldchange
		DESeq_d[ID]['log2change']= log2change
		DESeq_d[ID]['pval']= pval
		DESeq_d[ID]['padj']= padj
		
	    else: #calculate the value using the scaled pcs
		#first check if only pcs are listed for this gene:
		ID= fields[1]
		#check if any of the counts >1. If not, don't add
		#could also change this to a higher cutoff (>2) to be more stringent if needed
		include=False
		for i in repkeys:    
		    if countd[ID][i]>includeco:
			include=True
		
		if include==True:
		    A_r1= pc_normedd[ID][repkeys[0]]['normed']
		    A_r2= pc_normedd[ID][repkeys[1]]['normed']
		    B_r1= pc_normedd[ID][repkeys[2]]['normed']
		    B_r2= pc_normedd[ID][repkeys[3]]['normed']
		   
		    baseMeanA= np.mean([A_r1, A_r2])
		    baseMeanB= np.mean([B_r1, B_r2])
		    
		    foldchange= baseMeanB/baseMeanA
		    log2change= math.log(foldchange, 2)
		    
		    DESeq_d[ID]={}
		    DESeq_d[ID]['baseMeanA']= baseMeanA
		    DESeq_d[ID]['baseMeanB']= baseMeanB
		    DESeq_d[ID]['foldchange']= foldchange
		    DESeq_d[ID]['log2change']= log2change
		    DESeq_d[ID]['pval']= float(fields[7])
		    DESeq_d[ID]['padj']= float(fields[8])
    return DESeq_d

def getVals(d, gene, keys):
    if gene in d:
	vals=[d[gene][i] for i in keys]
    else:
	vals=['ND']*len(keys)
    return vals


def getTEvals(FPvals, Tvals):
    
    #get delta TE
    TEchange= FPvals[0]-Tvals[0]
    
    #get the baseValues, don't log transform for now, can do later.
    WT_TE= FPvals[3]/Tvals[3]
    mut_TE= FPvals[4]/Tvals[4]
    
    return [TEchange, WT_TE, mut_TE]

def parseCounts(countfile):
    '''Parse the count file, i.e. DESeq input.'''
    countdict={}
    
    f=open(countfile, 'r')
    libs= f.readline().strip('\n').split('\t')[1:]
    for i in libs:
        countdict[str(i)]={}
        
    for line in f:
        fields= line.strip('\n').split('\t')
        gene= fields[0]
        vals= fields[1:]
        countdict[gene]={}
        for i in range(0, len(vals)):
            countdict[gene][libs[i]]= float(vals[i])
    f.close()
    return countdict

def present(countd, gene):
    '''return true if there are values for at least one of the libraries. False otherwise'''
    for gene in countd:
        #determine if 0 in all of the libraries, if so, don't write to file
        hasgene=False
        for i in countd[gene]:
            if countd[gene][i]>0:
                hasgene=True
                
    return hasgene

def normCounts(countdict, sfactors):
    '''Given the size factors and the dictionary of counts, return normalized counts to use for visualization.
    Add a scaled pseudocount if there are zeros in some, but not all of the libraries. If no counts in any of the libraries, remove the gene'''
    
    normedd={}
    pc_normedd={}
    
    for gene in countdict:
	if present(countdict, gene): #Currently this checks to see if any of the libraries have any counts and then adds pcs if yes.
	    normedd[gene]={}
	    pc_normedd[gene]={}
	    for lib in countdict[gene]:
		normedd[gene][lib]={}
		pc_normedd[gene][lib]={}
		val= countdict[gene][lib]
		if val<1: #add a pseudocount for genes here
		    pc_val=1
		    npc_val= pc_val/sfactors[lib]
		    normedd[gene][lib]['normed']= val #these should both be 0 now
		    normedd[gene][lib]['counts']= val
		    pc_normedd[gene][lib]['normed']= npc_val
		    pc_normedd[gene][lib]['counts']= pc_val
		else:
		    normedval= val/sfactors[lib]
		    normedd[gene][lib]['normed']=normedval
		    normedd[gene][lib]['counts']= val
		    pc_normedd[gene][lib]['normed']= normedval #will be the same for genes not affected by the pseudocount
		    pc_normedd[gene][lib]['counts']= val
		
    return normedd, pc_normedd


def updateRpkmPickle(rpkmspickle, normedd):
    
    f= open(rpkmspickle, 'r')
    d= pickle.load(f)
    f.close()
    
    #add the normalized counts to the dictionary, also add '' where any missing values occur
    for gene in d:
	for lib in d[gene]:
	    if (gene in normedd) and (lib in normedd[gene]):
		d[gene][lib]['scaledcounts']= normedd[gene][lib]['normed'] #b/c I did not add the mt genome ones to DESeq, those will not have a 'scaledcounts' value
    
    f= open(rpkmspickle, 'w')
    pickle.dump(d, f)
    f.close()
	    
def writeOutPickle(rundict, originalD, countdict, normeddict, pc_normeddict, controlName='WT'):
    '''For the controlName, will usually be WT. Use the normeddict to see if we want to add these values.
    Use the pc_normeddict to see which values to add'''
    
    FPFile= rundict['FP']
    TFile= rundict['tot']
    outname= rundict['outname']
    basename= rundict['basename']
    
    Brep1Fkey=originalD[basename]['F'][0]
    Brep2Fkey=originalD[basename]['F'][1]
    Brep1Tkey=originalD[basename]['T'][0]
    Brep2Tkey=originalD[basename]['T'][1]
    
    Arep1Fkey=originalD[controlName]['F'][0]
    Arep2Fkey=originalD[controlName]['F'][1]
    Arep1Tkey=originalD[controlName]['T'][0]
    Arep2Tkey=originalD[controlName]['T'][1]
    
    FPrepkeys= [Arep1Fkey, Arep2Fkey, Brep1Fkey, Brep2Fkey] #use this to fetch the original values to calculate a value for the pseudocount cases
    Trepkeys= [Arep1Tkey, Arep2Tkey, Brep1Tkey, Brep2Tkey]
    
    FPd= parseDESeq(FPFile, countdict, pc_normeddict, FPrepkeys) #get vals from DE-Seq file
    Td= parseDESeq(TFile, countdict, pc_normeddict, Trepkeys)
    
    allgenes=set(FPd.keys()).union(set(Td.keys())) #some of the FP genes will not have matching total genes
        
    k=open(outname, 'w')
  
    d={}
    for gene in allgenes: #at this point they should already be filtered for ones occurring in at least one library
        d[gene]={}
        DEkeys=['log2change', 'pval', 'padj', 'baseMeanA', 'baseMeanB']
        TEkeys=['log2change']
        
        FPvals= getVals(FPd, gene, DEkeys)
        Tvals= getVals(Td, gene, DEkeys)
	
	if FPvals[0]!='ND':        
	    WT_fp1= pc_normeddict[gene][Arep1Fkey]['normed']
	    WT_fp2= pc_normeddict[gene][Arep2Fkey]['normed']
	    Mut_fp1= pc_normeddict[gene][Brep1Fkey]['normed']
	    Mut_fp2= pc_normeddict[gene][Brep2Fkey]['normed']
	    
	    WT_fp1_counts= normeddict[gene][Arep1Fkey]['counts'] #include the original, not pc values for each
	    WT_fp2_counts= normeddict[gene][Arep2Fkey]['counts']
	    Mut_fp1_counts= normeddict[gene][Brep1Fkey]['counts']
	    Mut_fp2_counts= normeddict[gene][Brep2Fkey]['counts']
	
	if Tvals[0]!='ND':
	
	    WT_t1= pc_normeddict[gene][Arep1Tkey]['normed']
	    WT_t2= pc_normeddict[gene][Arep2Tkey]['normed']
	    Mut_t1= pc_normeddict[gene][Brep1Tkey]['normed']
	    Mut_t2= pc_normeddict[gene][Brep2Tkey]['normed']
	    
	    WT_t1_counts= normeddict[gene][Arep1Tkey]['counts']
	    WT_t2_counts= normeddict[gene][Arep2Tkey]['counts']
	    Mut_t1_counts= normeddict[gene][Brep1Tkey]['counts']
	    Mut_t2_counts= normeddict[gene][Brep2Tkey]['counts']
	
	if (FPvals[0]!='ND') and (Tvals[0]!='ND'):
	
	    TEvals= getTEvals(FPvals, Tvals)
	    
	    d[gene]['TE']={}
	    d[gene]['TE']['val']= TEvals[0]
	    d[gene]['TE']['baseMeanA']= TEvals[1]
	    d[gene]['TE']['baseMeanB']= TEvals[2]
	    
	    #get the TEs from individual replicates, these shouldn't throw errors, b/c we already got rid of zero val ones
	    d[gene]['TE']['A_r1']= WT_fp1/WT_t1
	    d[gene]['TE']['A_r2']= WT_fp2/WT_t2
	    d[gene]['TE']['B_r1']= Mut_fp1/Mut_t1
	    d[gene]['TE']['B_r2']= Mut_fp2/Mut_t2
	
	if FPvals[0] !='ND':
	    d[gene]['FP']={}
	    d[gene]['FP']['val']= FPvals[0]
	    d[gene]['FP']['pval']= FPvals[1]
	    d[gene]['FP']['padj']= FPvals[2]
	    d[gene]['FP']['baseMeanA']= FPvals[3]
	    d[gene]['FP']['baseMeanB']= FPvals[4]
	    d[gene]['FP']['A_r1']= WT_fp1
	    d[gene]['FP']['A_r2']= WT_fp2
	    d[gene]['FP']['B_r1']= Mut_fp1
	    d[gene]['FP']['B_r2']= Mut_fp2
	    
	    d[gene]['FP']['A_counts_r1']= WT_fp1_counts
	    d[gene]['FP']['A_counts_r2']= WT_fp2_counts
	    d[gene]['FP']['B_counts_r1']= Mut_fp1_counts
	    d[gene]['FP']['B_counts_r2']= Mut_fp2_counts
	
	if Tvals[0]!='ND':
	    d[gene]['tot']={}
	    d[gene]['tot']['val']= Tvals[0]
	    d[gene]['tot']['pval']= Tvals[1]
	    d[gene]['tot']['padj']= Tvals[2]
	    d[gene]['tot']['baseMeanA']= Tvals[3]
	    d[gene]['tot']['baseMeanB']= Tvals[4]
	    d[gene]['tot']['A_r1']= WT_t1
	    d[gene]['tot']['A_r2']= WT_t2
	    d[gene]['tot']['B_r1']= Mut_t1
	    d[gene]['tot']['B_r2']= Mut_t2
	    
	    d[gene]['tot']['A_counts_r1']= WT_t1_counts
	    d[gene]['tot']['A_counts_r2']= WT_t2_counts
	    d[gene]['tot']['B_counts_r1']= Mut_t1_counts
	    d[gene]['tot']['B_counts_r2']= Mut_t2_counts
	    
    pickle.dump(d, k)
    k.close()
        
def main(argList):
    
    params= argList[0]
    
    #get the names of all the files needed:
    prefix= os.path.join(params.outdir, params.bigoutname)
    countfile= prefix+'_counts.txt'
    sfactorfile= prefix+'_sizefactors.txt'
    
    countdict= parseCounts(countfile)
    sfactors= parseSfactors(sfactorfile)
    
    #get the normalized values for each gene
    normedd, pc_normedd= normCounts(countdict, sfactors)
    
    #make a list of the analyses to run:
    toRun=[]
    originalD={}
    
    for exp in range(0, len(params.reporder)):
	WT= params.WT[exp]
	for strain in params.reporder[exp]:
	    originalD[strain]={}
	    originalD[strain]['F']= params.replicates[strain][0:2]
	    originalD[strain]['T']= params.replicates[strain][2:4]
	    if strain != WT:
		rund={}
		rund['basename']= strain
		rund['FP']= '%s_%s_Fanalysis.txt' % (prefix, strain)
		rund['tot']= '%s_%s_Tanalysis.txt' % (prefix, strain)
		rund['outname']= '%s_%s_allData.p' % (prefix, strain)
		rund['sdict']= sfactors
		rund['control']= WT
		
		toRun.append(rund)
	
    
    for i in toRun:
        writeOutPickle(i, originalD, countdict, normedd, pc_normedd, controlName= i['control']) #use pc to write the experiment output
    
    print 'updating the rpkm dict with the scaled values'
    #write the scaled data to the rpkms.p dict
    rpkm_pickle= '%s_rpkms.p' % os.path.join(params.outdir, params.bigoutname) #use the non-pc to write the the files for GEO
    updateRpkmPickle(rpkm_pickle, normedd)
    
    
if __name__ == '__main__':
    main(sys.argv[1])