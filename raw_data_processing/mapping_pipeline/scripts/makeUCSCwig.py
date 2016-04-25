#Boris Zinshteyn
#Makes a UCSC-compatible WIG file
#arguments:
#   pickledWeightedReads
#   annotationCounts - should have a line such that starts w/ "CDS", the  this value will be used to normalize all of the reads
#   genomePickle
#   outputPrefix
import sys, cPickle as pickle, os, gzip

def getCDSReadCount(annotationCounts):
    f=open(annotationCounts)
    for line in f:
        ll=line.strip().split('\t')
        if ll[0] == 'CDS':
            counts = float(ll[1])
    f.close()
    return counts


def main():


    pickledWeightedReads, annotationCounts, outputPrefix, name= sys.argv[1:]
    
    weightedReads = pickle.load(open(pickledWeightedReads))
    
    normalization = getCDSReadCount(annotationCounts)/1000000
    
    plusWig = gzip.open(outputPrefix+'_plus.wig.gz', 'w')
    minusWig = gzip.open(outputPrefix+'_minus.wig.gz', 'w')
    plusWig.write('track type=wiggle_0 name=%s\n' % (name+'_plus'))
    minusWig.write('track type=wiggle_0 name=%s\n' % (name+'_minus'))
    allChrs=set(weightedReads['+'].keys()).union(set(weightedReads['-'].keys()))
    for chr in allChrs:
        #minPos = min([min(weightedReads['+'][chr].keys()), min(weightedReads['-'][chr].keys())])
        #maxPos = max([max(weightedReads['+'][chr].keys()), max(weightedReads['-'][chr].keys())])
        plusWig.write('variableStep chrom=%s\n' % (chr))
        minusWig.write('variableStep chrom=%s\n' % (chr))
        if chr in weightedReads['+']:
            for i in sorted(weightedReads['+'][chr].keys()):
                plusWig.write('%d\t%f\n' % (i, weightedReads['+'][chr][i]/normalization))
        if chr in weightedReads['-']:
            for i in sorted(weightedReads['-'][chr].keys()):
                minusWig.write('%d\t%f\n' % (i, weightedReads['-'][chr][i]/normalization))
    plusWig.close()
    minusWig.close()
        
    
            
    


main()

#>>> print '%.2f' %answer
#0.83

