# ORF_length_comparisons

1) Parse other organism datasets to convert to common pickled dictionary format:
python parseGuo.py mouse_mm9_refseq.gtf GSM546987_fp_WT_quantification.txt GSM546988_mrna_WT_quantification.txt <outprefix> --writeLenP
python parseGuo.py human_hg18_refseq.gtf GSM546920_fp_mock_32hr_quantification.txt GSM546921_mrna_mock_32hr_quantification.txt <outprefix> --writeLenP
python parseStadler.py c_elegans.WS230.annotations.gff3 GSE48140_Celegans_RPF_deseq.txt GSE48140_Celegans_mRNA_deseq.txt <outprefix> --writeLenP


2) Make plots of TE vs. ORF length in other organisms:
python plotLenvTE_conservation.py Guo_mouse.p <outprefix> mouse
python plotLenvTE_conservation.py Guo_hela.p <outprefix> human
python plotLenvTE_conservation.py Stadler_worm_A.p <outprefix> worm

3) Look at mRNA attributes vs. deltaTE in M1X
python collectAtts.py <mRNAattributes.csv> Pelechano_2014_S3.txt Pelechano_5p_UTRs_MFEs.txt Pelechano_3p_UTRs_MFEs.txt  sacCer.gff <allData_filtered_file> <allData_file> <rpkm_file> <outprefix>
python plotAtts.py <mergedAttributes> <outprefix>
python parsedeGodoy.py deGodoy_2008_S5.csv <allData_filtered_file>

4) Evaluate ORF vs. transcript length for determining TE:
python partialCorr_wUTRs.py <mergedAttributes> <yeast_ORF_lens> <outprefix>

notes:
parseStadler.py produces two files, _A (diapause) and _B (fed), _A is presented in the supplemental figure

file descriptions:

<outname> prefix for output file or graph
<mRNAattributes.csv> is a csv file containing mRNA attributes provided by Eckhard Jankowsky,
only attritbute used from this file is tAI (tRNA adaptation index), which they calculated using values from Tuller et al., 2010

Pelechano_2014_S3.txt
Table S3 from Pelechano et al., 2014

Pelechano_5p_UTRs_MFEs.txt/ Pelechano_3p_UTRs_MFEs.txt
minimum free energy (MFE) of UTR sequences using median UTR lengths from Pelechano et al., 2014. Calculations done using RNAfold (Gruber et al., 2008)

sacCer.gff
Saccharomyces cerevisiae gff annotation file downloaded from Saccharomyces Genome Database Sept. 2, 2011

<allData_filtered_file> pickled python dictionary containing experiment data (e.g. yeast_mutants_allData_M1X_filtered.p)

<allData_file> pickled python dictionary containing experiment data, not filtered by read count (e.g. yeast_mutants_allData_M1X.p)

<rpkm_file> pickled python dictionary containing rpkms (i.e. yeast_mutants_rpkms.p)

<mergedAttributes> csv file output of collectAtts.py
e.g. provided file, 141104_filtered_final.csv

deGodoy_2008_S5.csv
Table S5 from de Godoy et al., 2008

<yeast_ORF_lens> pickled python dictionary mapping gene->['average', 'median', 'ORFLen', 'sepLen'], e.g. ORF_length_line/ORFlenpickles/150507_yeastLens.p
ORFLen= ORF lengths; average, median= average, median transcript lengths; sepLen= transcript length based on medians of 5' and 3' UTR separately.

other organism's annotation files:

downloaded from UCSC May 31, 2015:
mouse_mm9_refseq.gtf
human_hg18_refseq.gtf

downloaded from WormBase June 15, 2015:
c_elegans.WS230.annotations.gff3


