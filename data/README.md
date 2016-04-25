# data

gene_lists/
    contains files containing genelists used in plots, i.e. RP and MRP lists, closed loop and strong closed loop groups
    
allData_files.tar.gz
    compressed filtered and unfiltered pickled python dictionaries (~80 MB uncompressed) containing experiment data. Referenced as <all_Data_file> and <all_Data_filtered_file> in other README files
    
    descriptions:
    <all_Data_file> pickled python dictionary containing experiment data (e.g. yeast_mutants_allData_M1X.p).
    <all_Data_filtered_file> pickled python dictionary containing experiment data, filtered for the 128 counts cutoff criteria (e.g. yeast_mutants_allData_M1X_filtered.p).
