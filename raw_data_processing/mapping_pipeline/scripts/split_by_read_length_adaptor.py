#Author: MKT
#Author: MKT
#purpose: to parse bowtie output files (after they have been weighted based on multiple mappings with weight_repeats.py and convert to pickled python dictionary. This one hashes by the start position of the read so it will be fast to classify them into exons, introns, etc. with the annotation file.
#version history:
#created 9/?/09
#updated 1/26/10- changed it to hash by start with respect to strand and not the middle of the read. For - strand reads this means adding 34 to the start w.r.t the + strand for a 35 bp read.
#updated 6/3/10 to take start not as +34, b/c i'm only mapping the first 21 bp
#updated 6/18/10 to also count the length of the fragments
#6/22/10- changed - strand to start - i
#6/25/10- Boris: added plotting function
#6/27/10- Boris: now splits into several pickle and text files by length
#10/26/2011 - Boris: now takes read length as an argument and adjusts it
import matplotlib
matplotlib.use('Agg') #stupid hack to make it not crash since server has no visualization backend
import sys, pickle, numpy, matplotlib.pyplot as plt, gzip

strands= ['+','-']
complements= {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}

def create_subdict(taileddict, bin):
    
    taileddict[bin]= {}
    for strand in strands:
        taileddict[bin][strand]= {}
        for chr in chrs:
            taileddict[bin][strand][chr]= {}

    return taileddict

def plot_histogram(frag_dist, outputPrefix):
    plt.rcParams['pdf.fonttype'] = 42
    plus=frag_dist['+']
    minus=frag_dist['-']
    x_labels=plus.keys()
    x_labels.sort()
    x_labels=x_labels#toss the first bin, which is ambiguous length
    bar_locations=numpy.arange(len(x_labels))
    width=.5
    plus_vals=[]
    minus_vals=[]
    plusTotal=0
    minusTotal=0
    for bin in x_labels:
        plusTotal+=plus[bin]
        minusTotal+=minus[bin]
    for bin in x_labels:
        plus_vals.append(plus[bin]/float(plusTotal+minusTotal))
        minus_vals.append(-1*minus[bin]/float(minusTotal+plusTotal))
    plus_plot=plt.bar(bar_locations[20:], plus_vals[20:], width, color='blue')
    minus_plot=plt.bar(bar_locations[20:], minus_vals[20:], width, color='red')
    plt.ylabel('Fraction of Reads')
    plt.xlabel('Read Length')
    plt.xlim(20, 37)
    plt.title(outputPrefix.split('/')[-2])
    x_labels[0]='?' #indicate that length zero are really ambiguous length
    plt.xticks(bar_locations[20:]+width/2., x_labels[20:])
    plt.savefig(outputPrefix+'.pdf', transparency=True, format='pdf')


def findlength(seq, adaptor):
    for i in range(len(seq)):
        if adaptor.startswith(seq[i:]):
            return i
    return len(seq)

def main():
    mappedFile, mappingLength, readLength, outFile, outFolder= sys.argv[1:]
    outFilesbyLength={}# a mapping of read lengths to the appropriate file to save to
#create dictionary to store reads of the form, sorted_dict= {'+':{'chr1':{120:3, 150:2,....
    tailed_dict= {}
    mappingLength=int(mappingLength)
    readLength=int(readLength)
    f= gzip.open(mappedFile, 'r')
    
    frag_dist= {}

    frag_range= [0]
    for i in range(0, readLength+1):
        frag_range.append(i)
    for strand in strands:
        frag_dist[strand]= {}
        for i in frag_range:
            frag_dist[strand][i]= 0   
    i= 0
    unknown_chrs=0
    for line in f:
        i+=1
        fields= line.strip().split('\t')
        counts= float(fields[0].split('&')[1])
        chr= fields[2].strip()
        strand= fields[1]
        seq= fields[4]
        cigarString = fields[6]
        mappingLength = int(fields[7])
        if strand== '+':
            start= int(fields[3])
        else: 
            start= int(fields[3])+(mappingLength-1)

    
        frag_length= mappingLength
        frag_dist[strand][frag_length]+= counts
        if frag_length not in outFilesbyLength:
            outFilesbyLength[frag_length]=open(outFolder+'length'+str(frag_length)+'.wt', 'w')
        outFilesbyLength[frag_length].write(line)
            
    for strand in frag_dist:
        strand_tot= 0
        for bin in frag_range:
            strand_tot+= frag_dist[strand][bin]

        #for bin in frag_range:
        #    print 'percent in', str(bin), str(float(frag_dist[strand][bin])*100/strand_tot)
    f.close()
    plot_histogram(frag_dist, outFile)

main()
