# polysome_quantification

Scripts for quantifying polysome regions and plotting polysome profiles

1) Quantify P/M and 40S/60S ratios:
python quantifyPsomes.py <outname> <indir> <paramfile>

e.g.
python quantifyPsomes.py testquant example_files/ example_files/WT_params.txt 


2) Plot polysome profiles:
python plotGradient.py <outname> <label_file> <csv_file1> <csv_file2> <csv_file3> <...> -xdim 4 -ydim 4 --noAxes -maxA 1.0 -clip_left 2 -clip_right 70

e.g.
python plotGradient.py testgraph example_files/WT_labels.txt example_files/sample1.csv example_files/sample2.csv example_files/sample3.csv -xdim 4 -ydim 4 --noAxes -maxA 1.0 -clip_left 2 -clip_right 70

file descriptions:

<outname> prefix for output file or graph
<indir> directory containing csv files for quantification
<paramfile> file containing parameters 
<label_file> file containing sample names and other annotations for plot
<csv_file*> csv file output from Gradient Master