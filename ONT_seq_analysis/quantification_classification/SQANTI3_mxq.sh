##!/usr/bin/env python
#! /usr/bin/bash


# OUTPUT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC
OUTPUT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter
WORKING=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Scripts/

mkdir -p $OUTPUT
mxqsub -w $WORKING -o $OUTPUT/SQANTI3_filter.log -j 20 -m 100G -t 10h \
sh $WORKING/SQANTI3_filter_run.sh            




   
