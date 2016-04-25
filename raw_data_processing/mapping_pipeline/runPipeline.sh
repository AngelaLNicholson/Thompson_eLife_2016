#!/bin/bash

#Author: MKT

#usage: bash runpipeline.sh prefsFile.txt
#purpose: a shell script to run bowtie and all the python scripts for ribosome profiling.
#Input is a preferences file that contains variable names and their values in the form: VAR_NAME=10

#check to see if there are four arguments. else print usage syntax
if [ $# -ne 1 ]
then 
    echo "usage: "
    echo "    bash runPipeline.sh <preferences file>"
    exit 1
fi

#########################################################
#load preferences from file and set up names and genomes#
#########################################################


PARAMETERS=$1
bash $PARAMETERS
source $PARAMETERS
mkdir -p $OUTPUT_DIRECTORY
mkdir -p "${OUTPUT_DIRECTORY}output"
mkdir -p "${OUTPUT_DIRECTORY}qc"
mkdir -p "${OUTPUT_DIRECTORY}intermediates"
mkdir -p "${OUTPUT_DIRECTORY}rRNA"
name="$OUTPUT_DIRECTORY${OUTPUT_PREFIX}"
LOG_FILE="$OUTPUT_DIRECTORY${OUTPUT_PREFIX}.log"
rRNA_LOG="$OUTPUT_DIRECTORY${OUTPUT_PREFIX}_rRNA.log"
INTER_PREFIX="${OUTPUT_DIRECTORY}intermediates/${OUTPUT_PREFIX}" #for intermediate files
cp $PARAMETERS ${name}.parameters

if $USE_PREDEFINED_STRAIN ; then
    if [ "$PREDEFINED_STRAIN" == "s288c" ]; then
        #set up variables for s288c
        bowtie_index_genome="${SCRIPTS_PATH}bowtieindexes/20110902_sacCer3_bt2/20110902_sacCer3"
        bowtie_index_rRNA="${SCRIPTS_PATH}bowtieindexes/20110902_sacCer3_bt2/20110902_sacCer3_rRNA"
        #bowtie_index_splice_junctions="${SCRIPTS_PATH}bowtieindexes/S288C100729/S288Csplice"
        splice_junc_annotations="genome_features/20110902_sacCer3/20110902_sacCer3_verified_orfs.juncs"
        pickled_genome="genome_features/20110902_sacCer3/20110902_sacCer3.p"
    elif [ "$PREDEFINED_STRAIN" == "sigma" ]; then
        #set up variables for sigma
        bowtie_index_genome="${SCRIPTS_PATH}bowtieindexes/20110113_sigma_bt2/20110113_sigma"
        bowtie_index_rRNA="${SCRIPTS_PATH}bowtieindexes/sigma_1_30_2011/rRNA_sigma"
        #bowtie_index_splice_junctions="${SCRIPTS_PATH}bowtieindexes/sigma_1_30_2011/sigma_splice100623"
        splice_junc_annotations="genome_features/1-30-2011_sigma_feature_gffs/sigma_annotations-all-verified.juncs"
        pickled_genome="genome_features/sigmav7complete.p"
    
    elif [ "$PREDEFINED_STRAIN" == "sigma_SGD_2014" ]; then
        #set up variables for sigma
        bowtie_index_genome="${SCRIPTS_PATH}bowtieindexes/20140629_sigmaSGD/Sigma1278b_full"
        bowtie_index_rRNA="${SCRIPTS_PATH}bowtieindexes/20140629_sigmaSGD/rRNA_SGD"
        #bowtie_index_splice_junctions="${SCRIPTS_PATH}bowtieindexes/sigma_1_30_2011/sigma_splice100623"
        splice_junc_annotations="genome_features/20140629_SigmaSGD/Sigma_pipeline_files/Sigma1278b_full.juncs"
        pickled_genome="genome_features/20140629_SigmaSGD/Sigma_pipeline_files/Sigma1278b_full.p"
    
    elif [ "$PREDEFINED_STRAIN" == "Celegans" ]; then
        #set up variables for C. elegans
        bowtie_index_genome="${SCRIPTS_PATH}bowtieindexes/Celegans_WS242"
        bowtie_index_rRNA="${SCRIPTS_PATH}bowtieindexes/Celegans_rRNA"
        splice_junc_annotations="genome_features/20140819_Celegans/Celegans_WS242.juncs"
        pickled_genome="genome_features/20140819_Celegans/Celegans_WS242_genome.p"
        
    else
        echo "$PREDEFINED_STRAIN is not a valid strain, currently s288c and sigma are supported" >> $LOG_FILE
        exit 1
    fi
fi



echo [`date`]: "$RAW_SEQUENCES"
echo [`date`]: "$RAW_SEQUENCES">$LOG_FILE

echo [`date`]: "####Pipeline Version Info####">>$LOG_FILE
echo [`date`]: `svn info`>>$LOG_FILE
echo [`date`]: "####End Pipeline Version Info####">>$LOG_FILE
#########################
#  Mapping Reads        #
#########################
if $MAP_READS
then

    echo [`date`]: "Trimming Adaptors or poly-A">>$LOG_FILE
    #RAW_SEQUENCES should now be a zipped .gz file
    cutadapt --adapter $ADAPTOR_SEQUENCE --overlap 3 --minimum-length $MIN_MAPPING_LENGTH "$RAW_SEQUENCES" --output "${INTER_PREFIX}_trimmed.gz">>$LOG_FILE 2>>$LOG_FILE
    echo [`date`]: "Finished Trimming Adaptors or poly-A">>$LOG_FILE
    echo>>$LOG_FILE

    echo [`date`]: "Mapping to genome and defined splice junctions with TopHat">>$LOG_FILE
    tophat2 --output-dir "${OUTPUT_DIRECTORY}intermediates" --library-type fr-unstranded --no-novel-juncs --no-novel-indels --raw-juncs $splice_junc_annotations --num-threads 8 $bowtie_index_genome "${INTER_PREFIX}_trimmed.gz" 1>>$LOG_FILE 2>>$LOG_FILE
    echo [`date`]: "Finished mapping to genome and defined splice junctions with TopHat" >> $LOG_FILE
    
    echo [`date`]: "untrim reads, weigh repeats, convert  to bowtie1 mapping format">>$LOG_FILE
    if $UNIQUE_MAPPING_ONLY; then
        MAPQ_CUTOFF=50
    else
        MAPQ_CUTOFF=1
    fi
    samtools view -q $MAPQ_CUTOFF -h "${OUTPUT_DIRECTORY}intermediates/accepted_hits.bam" | python ${SCRIPTS_PATH}restoreReadsWeighRepeats.py ${INTER_PREFIX}_full_weighted_mappings.wt.gz 1>>$LOG_FILE 2>>$LOG_FILE
    echo [`date`]: "Finished untrim reads, weigh repeats, convert  to bowtie1 mapping format">>$LOG_FILE
fi

################################################
#  Analyze rRNA contamination (optional)      #
###############################################

if $ANALYZE_rRNA
then
    
    echo [`date`]: "Mapping reads to rRNA with bowtie1">>$LOG_FILE
    gunzip -c "${INTER_PREFIX}_trimmed.gz"| bowtie1 -k 6 -p 8 -v 3 $bowtie_index_rRNA - "${OUTPUT_DIRECTORY}rRNA/${OUTPUT_PREFIX}_rRNA.mp" 1>>$rRNA_LOG 2>>$rRNA_LOG
    
    echo [`date`]: "Plotting length and positions of rRNA fragments">>$LOG_FILE
    python ${SCRIPTS_PATH}analyzerRNA.py "${OUTPUT_DIRECTORY}rRNA/${OUTPUT_PREFIX}_rRNA.mp" "${OUTPUT_DIRECTORY}rRNA/" -species $species -oligos $oligos 1>>$rRNA_LOG 2>>$rRNA_LOG
    rm "${OUTPUT_DIRECTORY}rRNA/${OUTPUT_PREFIX}_rRNA.mp"
fi

###############################
# Launch Analyses             #
###############################
if $ANALYZE_COMBINED; then
    echo [`date`]: "##########launching analysis subpipeline combined reads############" >>$LOG_FILE
    bash analysis_subpipeline.sh $PARAMETERS "${INTER_PREFIX}_full_weighted_mappings.wt.gz" "${OUTPUT_DIRECTORY}" 1>>$LOG_FILE 2>>$LOG_FILE
fi
wait

if $ANALYZE_BY_LENGTH; then
	mkdir -p "${OUTPUT_DIRECTORY}by_length/"
    echo [`date`]: ${RAW_SEQUENCES##*/}:getting length dist... >>$LOG_FILE
    python ${SCRIPTS_PATH}split_by_read_length_adaptor.py ${INTER_PREFIX}_full_weighted_mappings.wt.gz $pickled_genome "${OUTPUT_DIRECTORY}output/${OUTPUT_PREFIX}_length_dist" "${OUTPUT_DIRECTORY}by_length/" 1>>$LOG_FILE 2>>$LOG_FILE
    #NOTE: make sure there's nothing in the 'by_length/' folder that contains .wt unless you want it analyzed
    echo [`date`]: "####data analysis by length bin######" >>$LOG_FILE
    for weightFile in $(ls "${OUTPUT_DIRECTORY}by_length/" | grep .wt.gz );
    do
        xpath=${weightFile%/*} 
        xbase=${weightFile##*/}
        xfext=${xbase##*.}
        xpref=${xbase%.*.*}
        echo [`date`]: "############${xpref}#############" >>$LOG_FILE
        bash analysis_subpipeline.sh $PARAMETERS "${OUTPUT_DIRECTORY}by_length/$weightFile" "${OUTPUT_DIRECTORY}by_length/${xpref}/" 1>>$LOG_FILE 2>>$LOG_FILE
    done

fi
wait



###############################
# Clean up intermediate Files #
###############################
if $CLEAN_UP; then
    echo [`date`] ${RAW_SEQUENCES##*/}: deleting intermediate files to save space >>$LOG_FILE
    rm "${INTER_PREFIX}_full_weighted_mappings.wt.gz"
    rm "${INTER_PREFIX}_trimmed.gz"
    rm "${OUTPUT_DIRECTORY}intermediates/accepted_hits.bam"
fi

echo [`date`] ${RAW_SEQUENCES##*/}: finished >>$LOG_FILE
echo [`date`] ${RAW_SEQUENCES##*/}: finished
