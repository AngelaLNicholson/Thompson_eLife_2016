# ORF_length_line

1) Plot 4G vs. ORF length
python makeORFLenLinePlots.py <outname> <yeast_ORF_lens> -datafiles <Park_4G_HP_file> -legends eIF4G -colors bluishGreen -datatypeKey TE --exclude_ASC1
#here the code was changed to exclude genes encoding eIF4G1 isoforms rather than Asc1 

2) Plot M1X vs. ORF length
python makeORFLenLinePlots.py <outname> <yeast_ORF_lens> -datafiles <allData_filtered_file> -legends M1X -colors vermillion -datatypeKey TE --exclude_ASC1

3) Plot M1X vs. ORF length +/- RPs:
python makeORFLenLinePlots.py <outname> <yeast_ORF_lens> -datafiles <allData_filtered_file> -legends M1X noRPs -colors vermillion blue -datatypeKey TE --exclude_ASC1 --exclude_RPs

file descriptions:

<outname> prefix for output graph

<all_Data_filtered_file> pickled python dictionary containing experiment data (e.g. yeast_mutants_allData_M1X_filtered.p).

<yeast_ORF_lens> pickled python dictionary mapping gene->['average', 'median', 'ORFLen', 'sepLen'], e.g. ORFlenpickles/150507_yeastLens.p
ORFLen= ORF lengths; average, median= average, median transcript lengths; sepLen= transcript length based on medians of 5' and 3' UTR separately.

Transcript length data from:
Pelechano et al., 2014

ORF length data extracted from Saccharomyces cerevisiae annotations gff downloaded from Saccharomyces Genome Database Sept. 2, 2011

Generation of <yeast_ORF_lens>:
python getRNALengths.py <outname>

<Park_4G_HP_file> pickled python dictionary containing eIF4G depletion data from Park et al., 2011. The heavy polysome data was used in the figure.

Generation of <Park_4G_HP_file>:
python make4Gpickle.py Park_4G_HP.csv Park_4G_LP.csv <outname>

Park_4G_HP.csv and Park_4G_LP.csv are files from Park et al., 2011, Table S1 and Table S3, respectively.