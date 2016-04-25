# Thompson_Elife_2016

Contents:

closedloop/
    Scripts used to examine attributes of the closed loop mRNAs (ORF lengths and TE changes in mutants)

common/
    Scripts used to generate plots from processed data, including histograms, violin, and scatter plots

data/
    Contains gene lists used in some analyses as well as experiment data files, noted as <allData_file> or <allData_filtered_file> in script documentation

GO_enrichment/
    Scripts used to test for GO category enrichment of affected mRNAs in mutants
  
halflife/
    Scripts used to examine and plot the relationship between TE and genome-wide measurements of mRNA halflives
 
ORF_length_comparisons/
    Scripts for correlating TE with ORF and transcript length, correlating TE changes in the asc1-M1X mutant with mRNA attributes, and correlating TE with ORF length in other organisms
  
ORF_length_line/
    Scripts used to plot the change in TE of mRNAs binned by ORF length

polysome_quantification/
    Scripts used to quantify polysomal regions measured by sucrose gradient fractionation. Also includes scripts for plotting polysome gradients.

raw_data_processing/
    Scripts used to process raw RNA-Seq reads, quantify counts for specific genes, and calculate ratios between different experimental samples
    
Version info:

Raw data processing was done using Tophat (v2.0.8b), cutadapt (1.4.2), samtools (0.1.18), bowtie (0.10.0.2), Bash (4.2.25), Python (2.7.3), R (3.1.2), run on Ubuntu 12.04.5 LTS
Downstream analysis and processing was performed using the Anaconda scientific computing package (Python 2.7.8 |Anaconda 2.1.0 (x86_64)|[GCC 4.2.1 (Apple Inc. build 5577)] on darwin)

