#Author: MKT
#purpose: to parse bowtie output files (after they have been weighted based on multiple mappings with weight_repeats.py and convert to pickled python dictionary. This one hashes by the start position of the read so it will be fast to classify them into exons, introns, etc. with the annotation file.
#version history:
#created 9/?/09
#updated 1/26/10- changed it to hash by start with respect to strand and not the middle of the read. For - strand reads this means adding 34 to the start w.r.t the + strand for a 35 bp read.
#4/1/10	Pavan: commented out superfluous print commands
#5/26/10 Boris: adjusting antisense read position by 20 instead of 34, since the mapping only uses 21 bp, this needs to be changed if the mapping size is changed. If read sizes aren't constant, then OH BOY!!! Also, added M1killer chromosome
#6/24/10 Boris: Corrected for bowtie vs genome indexing
#6/23/2013 Boris: Removed arbitrary hard-coding chromosome names, and accounted for split reads (indels, splice junctions)
strands= ['+','-']

import sys, cPickle as pickle, gzip, re

def parseCigarString(cigarString, mappingLength):
    '''
    arguments:
        cigarString - a string of [#][A-Z][#][A-Z] etc... that described the alignment of a read to the genome
    returns:
        the genomic distance spanned by this read
        
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch
    H can only be present as the
    rst and/or last operation.
    S may only have H operations between them and the ends of the CIGAR string.
    For mRNA-to-genome alignment, an N operation represents an intron. For other types of
    alignments, the interpretation of N is not de
    ned.
    4 Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ. 
    '''
    
    numbers = re.split('[A-Z,=]', cigarString)[:-1]#not sure why, bu this seems to produce an extra blank entry thay I'm stripping off
    tags = re.split('[0-9]*', cigarString)[1:]
    
    genomeMulitpliers = {'M':1, 'I':0, 'D':1, 'N':1, '=':1, 'X':1}
    readMulitpliers = {'M':1, 'I':1, 'D':0, 'N':0, '=':1, 'X':1}
    
    readMappingSpan = 0
    genomeMappingSpan = 0
    
    for i in range(len(tags)):
        tag = tags[i]
        genomeMulitplier = genomeMulitpliers[tag]
        readMulitplier = readMulitpliers[tag]
        tagLength = int(numbers[i])
        readMappingSpan += readMulitplier*tagLength
        genomeMappingSpan += genomeMulitplier*tagLength
        
    assert readMappingSpan == mappingLength
    return readMappingSpan, genomeMappingSpan

sortedfile, pickled_sorted= sys.argv[1:]

#create dictionary to store reads of the form, sorted_dict= {'+':{'chr1':{120:3, 150:2,....
sorted_dict= {}
for strand in strands:
    sorted_dict[strand]={}


f= gzip.open(sortedfile, 'r')
i= 0
for line in f:
    i+=1
    fields = line.strip().split('\t')
    counts = float(fields[0].split('&')[1])
    chr = fields[2].strip()
    strand = fields[1]
    cigarString = fields[6].strip()
    mappingLength = int(fields[7].strip())
    
    readMappingSpan, genomeMappingSpan= parseCigarString(cigarString, mappingLength)
    
    if strand== '+':
        start= int(fields[3])
    else: 
        start= int(fields[3])+genomeMappingSpan-1 #When a read maps to the minus strand, bowtie returns the reverse complement, and indicates where this reverse mapped on the + strand. Thus the original 5' end of the read actually was 20nt downstream on the + strand, since we are using the first 21nt for mapping. If that changes, this will give wrong results.
    
    if chr not in sorted_dict[strand]:
        sorted_dict[strand][chr]={}
    if start in sorted_dict[strand][chr]: #check to see if start has already been hashed into dict
        sorted_dict[strand][chr][start]+= counts #just add the number of counts to that start position

    else:
        sorted_dict[strand][chr][start]=counts #just need to keep track of the number of counts from that position. 
    


f.close()

g= open(pickled_sorted, 'w')
pickle.dump(sorted_dict, g)
g.close()

