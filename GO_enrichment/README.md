# GO_enrichment

1) build GO tree-> yeast genes file
python buildGOtree.py go-basic.obo gene_association.sgd.gz <outname> -format SGD

2) test GO category enrichment for mRNAs with altered FP, total, total in yeast mutants
python GOstats.py <outname> <yeastGOtree> -files <all_Data_file> -min_genes 21 -max_genes 300 --mw --write_numgenes

3) make plot of GO category ORF length vs. deltaTE
python plotGObyLen.py <yeastGOtree> <yeast_ORF_lens> <all_Data_filtered_file> <outname> 21

file descriptions:

<outname> prefix for output

go-basic.obo
GO annotation file

gene_association.sgd.gz
Gene associations with GO categories

Source for go-basic.obo and gene_association.sgd.gz, downloaded Oct. 4, 2013
http://www.geneontology.org/GO.downloads.annotations.shtml

<yeastGOtree> pickled python dictionary mapping GO id->associated genes, output from buildGOtree.py

<all_Data_file> pickled python dictionary containing experiment data (e.g. yeast_mutants_allData_M1X.p).

<all_Data_filtered_file> pickled python dictionary containing experiment data, filtered for the 128 counts cutoff criteria (e.g. yeast_mutants_allData_M1X_filtered.p).

<yeast_ORF_lens> pickled python dictionary mapping gene->['average', 'median', 'ORFLen', 'sepLen'], e.g. ORF_length_line/ORFlenpickles/150507_yeastLens.p
ORFLen= ORF lengths; average, median= average, median transcript lengths; sepLen= transcript length based on medians of 5' and 3' UTR separately.

Transcript length data from:
Pelechano et al., 2014

ORF length data extracted from Saccharomyces cerevisiae annotations gff downloaded from Saccharomyces Genome Database Sept. 2, 2011
