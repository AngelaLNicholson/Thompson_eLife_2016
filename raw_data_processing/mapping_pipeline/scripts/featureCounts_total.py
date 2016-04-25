#in order to count number of reads in introns, exons and intergenic regions. takes in dictionary with start positions as keys and number of reads as values
#i.e. dict= {+:{chr1:{1096:1, 456:3.5...}}}
#uses bowtie notation, i.e. +/-, two_micron

#6/2/10	Boris- modified to do all annotations at once using SGD gff file with all annotations
#2/21/2011  Boris - Adjusted script to offset all features by 14 nt, to account for p site being 12nt down from 5' end of reads, and 2 nt of fudge factor.
import sys, pickle
from collections import defaultdict

#chrs= ['chr1', 'chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chrmt','two_micron', 'M1Killer', 'LA']

#chrMap={'I':1, 'II':2, 'III':3, 'IV':4, 'V':5, 'VI':6, 'VII':7, 'VIII':8, 'IX':9, 'X':10, 'XI':11, 'XII':12, 'XIII':13, 'XIV':14, 'XV':15, 'XVI':16, 'Mito':'mt'}

strands= ['+','-']

##these should match the feature names in the s288c whole-genome .gff file from SGD
features=['CDS', 'intron', 'rRNA', '5_UTR', '3_UTR', 'tRNA']

def get_anns(annotations):

    ann_dict= defaultdict(lambda: defaultdict(list))
            
    g= open(annotations, 'r')

    for line in g:
        if not line.startswith('#'):
            orf_dict= {} #for each orf, will be of the form {name:..., start:..., end:....} #making a list of dictionaries. Is this inefficient or should I be hashing by something else??

            fields= line.strip().split()
            strand=fields[6]
            feature=fields[2]
            if strand in strands and feature in features:
                chr=fields[0]
                source=fields[1]
                start=fields[3]
                end=fields[4]
                #id=fields[8].split(';')[0].split('=')[1]

                #orf_dict['id']= id
                orf_dict['feature']=feature
                orf_dict['start']= int(start)
                orf_dict['end']= int(end)
                orf_dict['length']=int(end)-int(start)+1
                ann_dict[strand][chr].append(orf_dict)
        if line.startswith('##FASTA') or line.startswith('>'):
            break
    return ann_dict

def score_orfs(pickledreads, ann_dict):

    tot_counts= 0

    f= open(pickledreads,'r')
    
    read_dict= pickle.load(f) #load the dictionary of reads

    #count the total number of reads in the file
    for strand in strands:
        for chr in read_dict[strand]:
            for i in read_dict[strand][chr].keys():
                tot_counts+=read_dict[strand][chr][i]
            #tot_counts+= sum(read_dict[strand][chr].values())
    
    shortest_annotations={}
    for strand in strands: #hash for each strand
        shortest_annotations[strand]= defaultdict(dict)#for each location a read maps, will keep the length of the shortest feature found there and it's type
    
    featureDict={}#will be a mapping of features to number of reads in that feature
    for strand in strands: #now go through each ORF annotation systematically and look up reads that map in that region.
        for chr in read_dict[strand]:
            if chr in ann_dict[strand]:
                for orf in ann_dict[strand][chr]: #first get all the dictionaries with start between start and stop of CDS
                    left_boundary= orf['start']-14
                    right_boundary= orf['end']-14
                    included= range(left_boundary, right_boundary+1) #create a range for every possible start point of the reads, assume exon boundaries included
                    
                    for i in included:
                    #see what reads map to the region, and add the annotation to each of them made it i-1 b/c bowtie is 0-based
                        if i-1 in read_dict[strand][chr]:
                            num_reads = read_dict[strand][chr][i-1]
                            if not shortest_annotations[strand][chr].has_key(i):
                                shortest_annotations[strand][chr][i]={}
                                shortest_annotations[strand][chr][i]['feature']=orf['feature']
                                shortest_annotations[strand][chr][i]['length']=orf['length']
                                shortest_annotations[strand][chr][i]['reads']=num_reads
                            elif shortest_annotations[strand][chr][i]['length']>orf['length']:
                                shortest_annotations[strand][chr][i]['feature']=orf['feature']
                                shortest_annotations[strand][chr][i]['length']=orf['length']
                                shortest_annotations[strand][chr][i]['reads']=num_reads
    for strand in strands:
        for chr in shortest_annotations[strand]:
            for region in shortest_annotations[strand][chr].values():
                if region['feature'] not in featureDict.keys():
                    featureDict[region['feature']]=region['reads']
                else:
                    featureDict[region['feature']]+=region['reads']
    for feature in features:
        if feature not in featureDict:
            featureDict[feature]=0
    f.close()
    return tot_counts, featureDict
    
def main():

    pickledreads, annotations, results= sys.argv[1:]

    ann_dict= get_anns(annotations)

    #print ann_dict
    tot_counts, featureDict= score_orfs(pickledreads, ann_dict)
    r= open(results, 'w')
    r.write('totCounts'+'\t'+str(tot_counts)+'\n')
    for feature, count in featureDict.iteritems():
        r. write('%s\t%f\t%.4f\n' % (feature, count, count/float(tot_counts)))
    r.close()
main()
