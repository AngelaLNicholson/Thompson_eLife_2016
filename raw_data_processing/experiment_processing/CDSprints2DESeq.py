#140611 Mary K. Thompson
#purpose: to convert the CDSprints files into 1) rpkms.p which contains both RPKM and counts per gene. 2) DE-Seq input counts file 3) DE-Seq Rscript

import sys
import os
import argparse
import cPickle as pickle

def unPickle(dictFile):
    f=open(dictFile)
    dict=pickle.load(f)
    f.close()
    return dict

def pickleoutPut(positions, outputFileName):
    '''
    Pickles the given dictionary
    '''
    g= open(outputFileName+".p", 'w')
    pickle.dump(positions, g)
    g.close()

def isInt(string):
    try:
        int(string)
        return True
    except:
        return False

def isFloat(string):
    try:
        float(string)
        return True
    except:
        return False

def countReads(footPrintMapping, offset, minmaplen):
    ''''
    Offset is then number of codons on the 5' side of the gene that we are excluding. The M of RPKM is only determined from the counted reads
    min_mapping_len= genes must have at least this many nt in the mapping length to be included
    '''
    minmapkb= minmaplen/1000.0
    
    countDict={}
    rpkmDict={}
    totalBodyReads=0
    for gene in footPrintMapping:
        mappingLengthInKb=(footPrintMapping[gene]['length']-3*offset-6)/1000.0 #kb of exon for this gene
        if mappingLengthInKb<minmapkb:
            continue
        mapping=footPrintMapping[gene]['hits']
        bodyCounts=0
        for i in range(-12+(offset*3), footPrintMapping[gene]['length']-18):
            if i in mapping:
                bodyCounts+=mapping[i]
        totalBodyReads+=bodyCounts
        countDict[gene]=bodyCounts
    
    millionsOfMappingReads=totalBodyReads/1000000.0
    
    #go back and add the total values and RPKMs to each gene
    for gene in countDict:
        #subtract 6 from the mapping length to account for starting with the start codon in the P-site and stoping with the stop codon in the A site
        bodyRpkm=(countDict[gene]/mappingLengthInKb)/(millionsOfMappingReads)
        rpkmDict[gene]=bodyRpkm
    return rpkmDict, countDict

def combineDicts(rpkmsDict, countsDict):
    '''change the hash from dataset->gene to gene->dataset'''
    
    d={}
    #get the union of all genes from all datasets here
    genes=set(countsDict[countsDict.keys()[0]].keys())
    for i in countsDict:
        genes= genes.union(set(countsDict[i].keys()))
        
    for gene in genes:
        d[gene]={}
        for fileName in countsDict:
            if gene in countsDict[fileName]:
                d[gene][fileName]= {}
                d[gene][fileName]['counts']= countsDict[fileName][gene]
                d[gene][fileName]['rpkm']= rpkmsDict[fileName][gene]

            else:
                d[gene][fileName]= {}
                d[gene][fileName]['counts']=0
                d[gene][fileName]['rpkm']=0

    return d

def writeCountsoneFile(d, params, countfilename):
    
    reps= params.replicates
    strain_order= params.strain_order
        
    g= open(countfilename, 'w')
    firstline='%s' % 'gene'

    order=[]
    
    for mutant in strain_order: 
        
        MutFP_1= reps[mutant][0]
        MutFP_2= reps[mutant][1]
        Muttot_1= reps[mutant][2]
        Muttot_2= reps[mutant][3]
        
        mut_order= [MutFP_1, MutFP_2, Muttot_1, Muttot_2]
        order.extend(mut_order)
            
    #write the header:            
    for j in order:
        firstline+='\t%s' % j
    firstline+='\n'
    g.write(firstline)
    
    #write the values gene by gene:
    for gene in d:
        if gene.startswith('Q'): #this will remove all the ones encoded in the mt genome
            continue
        else:
            g.write('%s' % gene)
            for dataset in order:
                g.write('\t%s' % d[gene][dataset]['counts'])
            g.write('\n')
    g.close()
    
def writeRscript(strain_order, countfilename, params):
    '''write the Rscript to be used for input into DESeq'''
    
    #generate required strings:
    basename= countfilename.split('_counts.txt')[0]
    sizefactorfile= basename+'_sizefactors.txt'
    dispersionfile= basename+'_dispersion.txt'
    rscriptfile= basename+'.R'
    
    mutantlist=[]
    FPlist=[]
    Totlist=[]
    for i in strain_order:
        FPname= '%s.F' % i
        Tname= '%s.T' % i
        FPlist.append('"'+ FPname+ '"')
        Totlist.append('"'+ Tname + '"')
        for j in range(2):
            mutantlist.append('"'+FPname+'"')
        for j in range(2):
            mutantlist.append('"'+Tname+'"')
    
    f= open(rscriptfile, 'w')
    
    #keep strings that need filling here
    FirstRcommands={0:'countTable = read.table(file = "%s", header = TRUE, row.names = 1 )' % countfilename, 1:'experimentDesign = data.frame(row.names = colnames(countTable), condition = c(%s), libType = c(rep("single-end",%s)))' % (','.join(mutantlist), len(mutantlist)), \
               2:'condition = experimentDesign$condition', 3:'library("DESeq")', 4:'cds = newCountDataSet(countTable, condition)', 5:'cds = estimateSizeFactors(cds)', 6:'write.table(sizeFactors(cds), "%s", quote=FALSE, sep=' % sizefactorfile, \
               7:'cds = estimateDispersions(cds, method="per-condition")'}
    
    for i in range(0, len(FirstRcommands)):
        line= FirstRcommands[i]
        if i==6:
            line+= r'"\t")'

        f.write(line+'\n')
    
    NextRcommands={'plotDispersion':'plotDispEsts(cds, name=%s)', 'writeTable':'write.table(fData(cds), "%s", quote=FALSE, sep=' % dispersionfile,\
                   'binomialTest':'res = nbinomTest( cds, %s)\nplotMA(res)\nhist(res$pval, breaks=100, main="")\nwrite.table(res, file="%s", quote=FALSE, sep=', 'devoff':'dev.off()'}
    
    #plot the dispersion:
    for i in FPlist:
        line= NextRcommands['plotDispersion'] % i
        f.write(line+'\n')
        
    for i in Totlist:
        line= NextRcommands['plotDispersion'] % i
        f.write(line+'\n')
        
    #wirte the dispersion table
    line= NextRcommands['writeTable']
    line+= r'"\t")'

    f.write(line+'\n\n')
    
    #run the binomial tests
    
    WTs=[0]
    for i in range(0, len(params.reporder)-1):
        nextI= WTs[i]+len(params.reporder[i])
        WTs.append(nextI)
    WTs.extend([len(params.strain_order)])
    
    #next go through and write all the pairs:
    for i in range(0, len(WTs)-1):
        start= WTs[i]+1
        end= WTs[i+1]
        for j in range(start, end):
            WTindex= i
            mutindex=j
  
            mutstring= ','.join([FPlist[WTs[i]], FPlist[j]])
            mutname= FPlist[j].strip('"').split('.')[0]
            outname= '%s_%s_%s' % (basename, mutname, 'Fanalysis.txt')
            line= NextRcommands['binomialTest'] % (mutstring, outname)
            line+= r'"\t")'
            f.write(line+'\n\n')
        
            mutstring= ','.join([Totlist[WTs[i]], Totlist[j]])
            mutname= Totlist[j].strip('"').split('.')[0]
            outname= '%s_%s_%s' % (basename, mutname, 'Tanalysis.txt')
            line= NextRcommands['binomialTest'] % (mutstring, outname)
            line+= r'"\t")'
            f.write(line+'\n\n')
        
    #dev off
    line= NextRcommands['devoff']
    f.write(line+'\n')

def main(argList):
        
    params= argList[0]
  
    #get the cds file names:
    cdsfiles= params.libraries.keys()
    outfilename= params.bigoutname
    outdir= params.outdir #need to add this all of the basenames
    offset= params.codons
    minmaplen= params.minmaplen
    
    fulloutname=os.path.join(outdir, outfilename)
    
    #count the reads
    rpkmDicts={}
    countDicts={}
    print 'counting genes'
    for name in cdsfiles:
        filepath= params.libraries[name]
        footPrintMapping=unPickle(filepath)#get the source data dictionary
        rpkmDicts[name], countDicts[name]=countReads(footPrintMapping, offset, minmaplen)#count reads in the features
    
    d= combineDicts(rpkmDicts, countDicts)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    
    f=open(fulloutname+'_rpkms.p', 'w')
    pickle.dump(d, f)
    f.close()
    
    print 'writing DESeq file'
    #make the DESeq count input file, exclude the mito genome genes (they all start with 'Q' as these are likely to be contaminants)
    countfilename= '%s_counts.txt' % fulloutname

    writeCountsoneFile(d, params, countfilename)
    
    #make the other DESeq input file:
    writeRscript(params.strain_order, countfilename, params)
    
if __name__ == '__main__':
    main(sys.argv[1:])