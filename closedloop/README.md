# closedloop

1) Make histogram of deltaTEs for closed loop mRNA groups in asc1-M1X mutant
python CostelloRegression.py <Costello_S5> <yeast_ORF_lens> <allData_filtered_file> <outname> --percentage -beta 5 --mw --plotlens

2) Make histogram of delta polysome association after eIF4G depletion
python CostelloRegression.py <Costello_S5> <yeast_ORF_lens> <allData_filtered_file> <outname> --percentage -beta 5 --mw --plotlens

3) Make histogram of ORF lengths among closed loop mRNAs
This plot is generated with either call above (1 or 2).

file descriptions:

<Costello_S5> Supplemental file 5 from Costello et al., 2015 in csv format

<outname> prefix for output graph

<yeast_ORF_lens> pickled python dictionary mapping gene->['average', 'median', 'ORFLen', 'sepLen'], e.g. ORF_length_line/ORFlenpickles/150507_yeastLens.p
ORFLen= ORF lengths; average, median= average, median transcript lengths; sepLen= transcript length based on medians of 5' and 3' UTR separately.

Transcript length data from:
Pelechano et al., 2014

ORF length data extracted from Saccharomyces cerevisiae annotations gff downloaded from Saccharomyces Genome Database Sept. 2, 2011

Generation of <yeast_ORF_lens>:
python getRNALengths.py <outname>
*see ORF_length_line/

<all_Data_filtered_file> pickled python dictionary containing experiment data (e.g. yeast_mutants_allData_M1X_filtered.p).

<Park_4G_HP_file> pickled python dictionary containing eIF4G depletion data from Park et al., 2011. The heavy polysome data was used in the figure.