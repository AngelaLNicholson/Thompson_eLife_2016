# experiment_processing

This folder contains scripts used to convert the counted reads into gene-specific counts, scale values using the DE-Seq package (hence the name runDESeq.py), and calculate ratios between experiment samples (in order to calculate TE values or changes in the experimental condition).

usage:

python runDESeq.py <parameter_file>

file descriptions:

<parameter_file> contains parameters for the analysis as well as locations for input and output on the machine, see example_input/runDESeq_example_params.txt