#!/bin/bash

#Author: Boris Zinshteyn
#purpose: just does the analysis part of the whole pipeline, created to allow pipeline to split reads by length and then analyze them
#and total RNA libraries

#version history:
#6/27/10- created
#1/30/2011 - added option for different yeast strains, currently s288c or sigma
#10/26/2011 - rewrote to use same paramters file as runpipeline.sh

#check to see if there are two arguments. else print usage syntax
if [ $# -ne 3 ]
then 
    echo "usage: "
    echo "     analysis_subpipeline <parameters file> </path/to/weighted_read_file.wt> <output dir>"
    echo 
    exit 1
fi


#take in arguments to the file: #1: scripts folder path and #2: raw dataset from the command line 3: strain name
PARAMETERS=$1
sh $PARAMETERS
source $PARAMETERS
weighted_reads=$2
OUTPUT_DIRECTORY=$3
mkdir -p $OUTPUT_DIRECTORY
mkdir -p ${OUTPUT_DIRECTORY}output
mkdir -p ${OUTPUT_DIRECTORY}qc
name="${weighted_reads%.*}"
outDir="${weighted_reads%/*}"
fileName="${name##*/}"
if $USE_PREDEFINED_STRAIN ; then
    if [ "$PREDEFINED_STRAIN" == "s288c" ]
    then 
    #set up variables for s288c
    annotations="genome_features/20110902_sacCer3/20110902_sacCer3.tsv"
    annotations_gff="genome_features/20110902_sacCer3/20110902_sacCer3.gff"
    cds_sequences="genome_features/20110902_sacCer3/20110902_sacCer3_coding_seqs.p"
    pickled_genome="genome_features/20110902_sacCer3/20110902_sacCer3.p"
    
    elif [ "$PREDEFINED_STRAIN" == "sigma" ]
    then
    #set up variables for sigma
    annotations="genome_features/tabbed verified ORF annotations/sigma_annotations-all-verified.tab"
    annotations_gff="genome_features/tabbed verified ORF annotations/sigmav7.gff"
    cds_sequences="genome_features/1-30-2011_sigma_feature_seqs/sigma_verified_coding.p"
    pickled_genome="genome_features/sigmav7complete.p"
    
    elif [ "$PREDEFINED_STRAIN" == "sigma_SGD_2014" ]
    then
    annotations="genome_features/20140629_SigmaSGD/Sigma_pipeline_files/Sigma1278b_full_verified_unchar.tab"
    annotations_gff="genome_features/20140629_SigmaSGD/Sigma_pipeline_files/Sigma1278b_full.gff"
    cds_sequences="genome_features/20140629_SigmaSGD/Sigma_pipeline_files/Sigma1278b_full_verified_unchar.p"
    pickled_genome="genome_features/20140629_SigmaSGD/Sigma_pipeline_files/Sigma1278b_full.p"
    
    elif [ "$PREDEFINED_STRAIN" == "Celegans" ]
    then
    annotations="genome_features/20140819_Celegans/Celegans_WS242.tab"
    annotations_gff="genome_features/20140819_Celegans/c_elegans.PRJNA13758.WS242.annotations.gff2"
    cds_sequences="genome_features/20140819_Celegans/Celegans_WS242_CDS.p"
    pickled_genome="genome_features/20140819_Celegans/Sigma_pipeline_files/Celegans_WS242_genome.p"
    
    else
        echo "$PREDEFINED_STRAIN is not a valid strain, currently s288c and sigma are supported"
        exit 1
    fi
fi
echo [`date`]: "pickling weighted reads"
#create database (in form of pickled python dictionary) to keep track of all the reads. Since Tophat output is 1-indexed now, no need to adjust indexing anymore
python ${SCRIPTS_PATH}pickleweighted.py "${weighted_reads}" "${OUTPUT_DIRECTORY}output/${fileName}.p"

echo [`date`]: "making cdsPrints"
# - Map footprints to each CDS, using coords relative to each start codon, returns a pickled dictionary
python ${SCRIPTS_PATH}mapfootprintstoCDS.py "${OUTPUT_DIRECTORY}output/${fileName}.p" "$annotations" "${OUTPUT_DIRECTORY}output/${fileName}_CDSprints" $FLANKING

###############################
#       Quality Control       #
###############################
echo [`date`]: "counting feature reads"
#count how many reads map to feature types defined in gff file. Type of features to count defined in script
python ${SCRIPTS_PATH}featureCounts_total.py "${OUTPUT_DIRECTORY}output/${fileName}.p" "$annotations_gff" "${OUTPUT_DIRECTORY}qc/${fileName}.annotation_counts" 

###############################
# Data Analysis               #
###############################
echo [`date`]: "plotting metaCDS"
python ${SCRIPTS_PATH}plotMetaCDS.py "${OUTPUT_DIRECTORY}output/${fileName}_CDSprints.p" $START_MIN $START_MAX $STOP_MIN $STOP_MAX "${OUTPUT_DIRECTORY}output/${fileName}_CDSplot" &

#make .wig file
echo [`date`]: "making .WIGs"
python ${SCRIPTS_PATH}makeUCSCwig.py "${OUTPUT_DIRECTORY}output/${fileName}.p" "${OUTPUT_DIRECTORY}qc/${fileName}.annotation_counts" "${OUTPUT_DIRECTORY}output/${fileName}" $fileName &

wait