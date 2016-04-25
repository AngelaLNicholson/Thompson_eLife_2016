#140630 Mary K. Thompson
#purpose: create directory with compressed files and md5sums to upload to GEO repository

import sys
import cPickle as pickle
import csv
import os
import gzip
import subprocess

def writeGEOfiles(d, sampleOrder, prefix):
    '''Write files suitable for upload to GEO'''
    
    f= gzip.open('%s_rpkms.txt.gz' % prefix, 'w')
    g= gzip.open('%s_counts.txt.gz' % prefix, 'w')
    h= gzip.open('%s_scaledcounts.txt.gz' % prefix, 'w')
    
    header='%s\t' % 'gene'
    header+= '\t'.join(sampleOrder)+'\n'
    f.write(header)
    g.write(header)
    h.write(header)
    
    #just round to one decimal place to save space
    for gene in d:
        #determine if 0 in all of the libraries, if so, don't write to file. These should be filtered out before here.
        hasgene=False
        for i in sampleOrder:
            if d[gene][i]['counts']>0:
                hasgene=True
                
        if hasgene==True:
            rpkmsline='%s\t' % gene
            countsline='%s\t' % gene
            rpkmsline+= '\t'.join(['%0.1f' % d[gene][i]['rpkm'] for i in sampleOrder])+'\n'
            countsline+= '\t'.join(['%0.1f' % d[gene][i]['counts'] for i in sampleOrder])+'\n'
            f.write(rpkmsline)
            g.write(countsline)
            
            #determine if the gene has scaledcounts. This will not be true for the mt genome genes, but each lib should be the same
            if 'scaledcounts' in d[gene][sampleOrder[0]]:
                scaledcountsline='%s\t' % gene
                scaledcountsline+='\t'.join(['%0.1f' % d[gene][i]['scaledcounts'] for i in sampleOrder])+'\n'
                h.write(scaledcountsline)
    f.close()
    g.close()
    h.close()

def main(argList):
    params= argList[0]
    
    outfilename= params.bigoutname
    outdir= params.outdir #need to add this all of the basenames
    
    geo_outdir= params.geo_outdir
    
    if not os.path.exists(geo_outdir):
        os.makedirs(geo_outdir)
   
    rpkm_pickle= '%s_rpkms.p' % os.path.join(outdir, outfilename)
    f= open(rpkm_pickle, 'r')
    d= pickle.load(f)
    f.close()
    
    sampleOrder= params.lib_order
    extraSamples= params.extra_order
    
    outprefix= os.path.join(geo_outdir, outfilename)
    writeGEOfiles(d, sampleOrder, outprefix)
    
    #gzip and make the GEO spreadsheet containing the md5sums
    geo_summary_file='%s_geosummary.txt' % outprefix
    geo_processed_md5_file='%s_processed_md5.txt' % outprefix
    geo_raw_md5_file='%s_raw_md5.txt' % outprefix
    
    summaryf= open(geo_summary_file, 'w')
    processedf= open(geo_processed_md5_file, 'w')
    rawf= open(geo_raw_md5_file, 'w')
    
    #the names of the processed files to write to the output    
    rpkmfile='%s_rpkms.txt.gz' % outprefix
    countsfile='%s_counts.txt.gz' % outprefix
    scaledcountsfile='%s_scaledcounts.txt.gz' % outprefix
            
    #going to hard link the raw files into the for geo submission folder, ln <oldfile> <newfile>
    allSamples=sampleOrder
    allSamples.extend(extraSamples)
    for i in allSamples:
        if i not in params.bmcnames: #this will prevent us from linking ones we haven't specified
            continue
        #get the raw file order:
        newRawName= os.path.join(geo_outdir, ('%s.fastq.gz' % i))
        oldRawName= os.path.join(params.bmcdir, params.bmcnames[i])
        if not os.path.exists(newRawName):
            ln_cmd= ['ln', oldRawName, newRawName] #link it if not already linked
            subprocess.check_call(ln_cmd)
        md5_cmd= ['md5sum', newRawName]
        result=subprocess.check_output(md5_cmd)
        md5result= result.split(' ')[0]
        
        rawf.write('%s\t%s\t%s\n' % (os.path.basename(newRawName), 'fastq', md5result))
        summaryf.write('%s\t%s\t%s\t%s\t%s\n' % (i, os.path.basename(newRawName), os.path.basename(rpkmfile), os.path.basename(countsfile), os.path.basename(scaledcountsfile)))
    
    #get the processed file md5sums
    pfiles=[rpkmfile, countsfile, scaledcountsfile]
    ptypes=['rpkms', 'raw counts', 'scaled counts']
    
    for i in range(0, len(pfiles)):
        md5_cmd=['md5sum', pfiles[i]] 
        result=subprocess.check_output(md5_cmd)
        md5result= result.split(' ')[0]
        processedf.write('%s\t%s\t%s\n' % (os.path.basename(pfiles[i]), ptypes[i], md5result))
    
    summaryf.close()
    processedf.close()
    rawf.close()
    
if __name__ == '__main__':
    main(sys.argv[1])
