# mapping pipeline

This folder contains scripts and sequence resources used to map the raw reads to the genome, count mapped reads and assign to features, store data in Python pickle (CDSprints.p files).

usage:

bash runPipeline.sh <parameter_file>

or

bash multiRunPipeline.sh <parameter_file1> <parameter_file2> <paramter_file3> <...>

file descriptions:

<parameter_file> contains parameters for the analysis as well as locations for input and output on the machine, see example prefs files in example_pref_files/