#! /usr/bin/python2.6
"""
mapfootprintstoCDS.py
Author: Boris Zinshteyn
purpose:To map the 5' ends of mapped reads relative to each ORF.
usage: python mapfootprintstoCDS.py pickled_reads.p tabbed_annotations outputPrefix flanking_nucleotides
input:
    pickledReads.p - a pickled read file as output by pickleweighted.py {strand:{chr:{position:#reads}}
    tabbbed_annotations - plain text, with gene    chr    strand    start    end    num_exons    exon_starts    exon_ends,   with stop always larger than start for a given exon (so - strand orientation needs to be flipped)
    outputPrefix - what you want the output called
    flanking_nucleotides - how many nucleotides on each side of each CDS to include in the resulting dictionary.
outputs:A pickled dictionary of {gene_name:{'chr':chr, 'strand':strand, 'start':gene_start, 'end':gene_end, 'length':feature_length, 'hits':{position_along_gene (a of AUG is 0):#_of_read_5'ends_mapping_here}

#Also counts how many reads are in each reading frame, and outputs that data to a pickled dictionary as well.
#i.e. dict= {+:{chr1:{1096:1, 456:3.5...}}}
#uses bowtie notation, i.e. +/-, two_micron
"""
import sys, pickle, numpy, math

strands= ['+','-']

"""
Input:
    annotations - the location of a tab-delimited text file containing: gene\tchr\tstrand\tstart\end\exon+starts\texon_ends

Output:
    ann_dict - a dictionary of {name:..., start:..., end:....}
"""
def get_anns(annotations):
##########################

    #make a dictionary to store the features in the annotation file ={'+':{'chr1':[orf1, orf2, orf3,...]..

    ann_dict= {}
    for strand in strands:
        ann_dict[strand]= {}
#        for chr in chrs:
#            ann_dict[strand][chr]= [] #make a list for each chromosome to store the ORFs
            
    g= open(annotations, 'r')

    for line in g:
        orf_dict= {} #for each orf, will be of the form {name:..., start:..., end:....} #making a list of dictionaries.

        fields= line.strip().split()
        id= fields[0]
        ids=id.split(',')

        if len(ids)> 1:
            print ids
            continue

        chr= fields[1]
        strand= fields[2]
        start= fields[3]
        end= fields[4]
        exons_start= []

        exons1= fields[6].split(',') #list of exon start coordinates as strings ['345', '789']
        for i in exons1:
            exons_start.append(int(i))

        exons_end= []

        exons2= fields[7].split(',') #list of exon end coordinates as strings
        for j in exons2:
            exons_end.append(int(j))
        
        orf_dict['id']= id
        orf_dict['exons_start']= exons_start
        orf_dict['exons_end']= exons_end
        if chr not in ann_dict[strand]:
            ann_dict[strand][chr]=[]
        ann_dict[strand][chr].append(orf_dict)
    return ann_dict


################################################
def map_cds_positions(pickledreads, ann_dict, libname, extension):
################################################
    
    ###########
    #Initialize dictionaries for counting total numbers of reads
    ###########
    num_footprints= 0
    cds_counts= 0
    #intron_counts= 0
    countedreads= 0
    
    ##########
    #Initialize dictionaries for storing analyzed data
    ##########
    counted_orfs= {} #make dictionary to keep track of all the reads mapping to each ORF
    frames={0:0, 1:0, 2:0}
    
    ######
    #Import source data from pickled dictionary
    #######
    f= open(pickledreads,'r')
    read_dict= pickle.load(f) #load the dictionary of reads
    
    ########
    #Go through each ORF annotation systematically and look up reads that map in that region,
    #and save them to a dictionary indexed with the first nucleotide of each gene's start codon as zero
    ########
    
    for strand in strands: 
        for chr in ann_dict[strand]:
            if chr in read_dict[strand]:
                for orf in ann_dict[strand][chr]: #first get all the dictionaries with start between start and stop of CDS
                    name= orf['id']
                    
                    #a mapping of read start sites in relations to the CDS start site
                    #(introns not included, so these are distances along the cDNA), to weighted number of reads starting at that location.
                    mapping_reads = {}
                    num_mapping_reads=0.0 #read numbers can be weighted, so need not be integers
                    num_exons= len(orf['exons_start'])
                    exon_start_offset=0 #the distance from the beginning of the mRNA to the start of this exon
                    
                    if strand=='+':
                        gene_start=orf['exons_start'][0]
                        gene_end=orf['exons_end'][num_exons-1]
                        for exon_index in range(num_exons):
                            exon_start=orf['exons_start'][exon_index]
                            exon_end=orf['exons_end'][exon_index]
                            exon_length=exon_end-exon_start+1
                            
                            left_extension, right_extension = 0, 0

                            #7/27/2011 BZ - should only have extensions for first and last exons, and only in proper directions
                            if exon_index==0:#first exon
                                left_extension += extension
                            if exon_index==num_exons-1:#last exon
                                right_extension += extension
                            
                            #loop through every position in the exon and any applicable extensions.
                            #if position appears in read dictionary, add the number of reads mapping there to the position
                            for i in range(exon_start-left_extension, exon_end+1+right_extension):
                                if i in read_dict[strand][chr]:
                                    rel_exon_pos=i-exon_start
                                    rel_cds_pos=rel_exon_pos+exon_start_offset
                                    if rel_cds_pos not in mapping_reads:
                                        mapping_reads[rel_cds_pos]=0
                                    mapping_reads[rel_cds_pos]+=read_dict[strand][chr][i]
                                    num_mapping_reads+=read_dict[strand][chr][i]
                                    frames[rel_cds_pos%3]+=read_dict[strand][chr][i]
                            exon_start_offset+=exon_length
                    elif strand=='-':
                        gene_start=orf['exons_end'][num_exons-1]
                        gene_end=orf['exons_start'][0]
                        for exon_index in range(num_exons-1, -1, -1):
                            exon_start=orf['exons_end'][exon_index]
                            exon_end=orf['exons_start'][exon_index]
                            exon_length=exon_start-exon_end+1
                            
                            #7/27/2011 BZ - should only have extensions for first and last exons, and only in proper directions
                            left_extension, right_extension = 0, 0
                            
                            if exon_index==0:#last exon
                                right_extension += extension
                            if exon_index==num_exons-1:#first exon
                                left_extension += extension
                                
                            #recall that in the annotations, genes on the - strand have the start and end indicated by the genomic coordinates, not transcript order,
                            #thus they need to be flipped, since we want to iterate through nucleotides 5'-->3' for clarity.
                            #so the exon starts will be larger numbers that exon ends, which is why we iterate backwards (by -1)
                            for i in range(exon_start+left_extension, exon_end-1-right_extension, -1):
                                if i in read_dict[strand][chr]:
                                    rel_exon_pos=exon_start-i
                                    rel_cds_pos=rel_exon_pos+exon_start_offset
                                    if rel_cds_pos not in mapping_reads:
                                        mapping_reads[rel_cds_pos]=0
                                    mapping_reads[rel_cds_pos]+=read_dict[strand][chr][i]
                                    num_mapping_reads+=read_dict[strand][chr][i]
                                    frames[rel_cds_pos%3]+=read_dict[strand][chr][i]
                            exon_start_offset+=exon_length
                    else:
                        print "unkown strand:"+strand
                    feature_length=exon_start_offset
                    counted_orfs[name]= {}
                    counted_orfs[name]['chr']=chr
                    counted_orfs[name]['strand']=strand
                    counted_orfs[name]['start']=gene_start
                    counted_orfs[name]['end']=gene_end
                    counted_orfs[name]['length']= feature_length
                    counted_orfs[name]['hits']= mapping_reads

    return counted_orfs, frames


def writeoutputfile(cds_positions, outputFileName):
    '''
    Prints the contents of the positions dictionary to to file in a tab-delimited format.
    '''
    try:
        outputFile=open(outputFileName, 'w')
    except:
        print "bad output"
        sys.exit()
    for orf in cds_positions.keys():
        hitsLine=''
        for position, reads  in cds_positions[orf]['hits'].iteritems():
            hitsLine+=(str(position)+':'+str(reads).strip()+',')
        hitsLine=hitsLine[:len(hitsLine)]
        outLine="%s\t%s\t%s\t%d\t%d\t%d\t%s\n" % (orf, cds_positions[orf]['chr'], cds_positions[orf]['strand'], cds_positions[orf]['start'], cds_positions[orf]['end'], cds_positions[orf]['length'], hitsLine)
        outputFile.write(outLine)
    outputFile.close()

def pickleoutPut(positions, outputFileName):
    '''
    Pickles the given dictionary
    '''
    g= open(outputFileName+".p", 'w')
    pickle.dump(positions, g)
    g.close()

######
#Appends the frame data from this analysis to a specified file
######
def appendFrame(frameFile, frames, outputFileName):
    outName=outputFileName.split('/')[2]+'_'+outputFileName.split('/')[-1].split('_')[-2]
    f=open(frameFile, 'a')
    outLine=outName
    total=float(sum(frames.values()))
    for frame in range(3):
        outLine+='\t'+str(frames[frame]/total)
    outLine+='\n'
    f.write(outLine)
    f.close
        
        
    
def main():

    reads, annotationfile, outputFileName, extension= sys.argv[1:]

    basename= reads.split('.')[0]
    
    ann_dict= get_anns(annotationfile)
#	num_footprints, intron_counts, exon_counts, inter_counts= score_orfs(pickledreads, ann_dict)

    cds_positions, frames = map_cds_positions(reads, ann_dict, basename, int(extension))
    
    #writeoutputfile(cds_positions, outputFileName)
    
    pickleoutPut(cds_positions, outputFileName)
    pickleoutPut(frames, outputFileName+'_frame')
    
    #appendFrame(frameFile, frames, outputFileName)

main()

