# raw_data_processing

Contents:

experiment_processing/
    Scripts used to convert the counted reads into gene-specific counts, scale values using the DE-Seq package (hence the name runDESeq.py), and calculate ratios between experiment samples (in order to calculate TE values or changes in the experimental condition).

mapping_pipeline/
    Scripts and sequence resources used to map the raw reads to the genome, count mapped reads and assign to features, store data in Python pickle (CDSprints.p files).

