##!/usr/bin/env python
#! /usr/bin/bash

DATA=/project/owlmayerTemporary/Derek/NGN3_project/RNA_alignGenome/alignment
WORKING=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Scripts

RESULT=/project/NeurodifferentiationBAF/Analyses_Derek/Results

mkdir -p $WORKING/mxq_log/AlignToOntTrans
 
mxqsub -w $WORKING -o $WORKING/mxq_log/AlignToOntTrans/salmon_nanopore.filtered_index.log -j 10 -m 50G -t 10h sh salmon_run.sh 


# for subdir in $DATA/day*/
# do
# 	# filename=$(basename -- "$f")
# 	# extension="${filename##*.}"
# 	# filename="${filename%.*}"

# 	subdir_name=$(basename "$subdir")
#     # echo "Subdirectory name: $subdir_name"

    
    
# 	mxqsub -w $WORKING -o $WORKING/mxq_log/"$subdir_name"_salmon.log -j 10 -m 50G -t 10h sh salmon_run.sh "$subdir_name"

# 	# echo $BAM/$filename.coverage
# 	# echo "$subdir_name"
# done