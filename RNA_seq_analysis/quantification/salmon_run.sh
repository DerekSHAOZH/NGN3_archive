##!/usr/bin/env python
#! /usr/bin/bash

source /home/shao/Tools/miniconda3/etc/profile.d/conda.sh
conda activate salmon

echo "$(date): Salmon starts"

## Indexing

# TRANSCRIPTOME=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_corrected.fasta
TRANSCRIPTOME=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Data/nanopore.filtered.fa

# RESULT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Data/nanopore_corrected_index
RESULT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Data/nanopore.filtered_index


salmon index -t $TRANSCRIPTOME -i $RESULT -p 10






# Quantification
# Check if a sample id is provided as an argument
# if [ -z "$1" ]; then
#     echo "Usage: $0 <sample_id>"
#     exit 1
# fi


# # Store the sample id from the command line argument
# SAMPLE=$1

# DATA=/project/owlmayerTemporary/Derek/NGN3_project/RNA_alignGenome/alignment

# TRANSCRIPTOME=/project/Neurodifferentiation_System/GeneralResources/gencode.v32.transcripts.fa.gz
# ANNOTATION=/project/Neurodifferentiation_System/GeneralResources/gencode.v32.primary_assembly.annotation.gtf
# INDEX=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Data/nanopore_corrected_index

# RESULT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans

# # echo $DATA/$SAMPLE
 

# # Alignment-based
# # salmon quant \
# #         -t $TRANSCRIPTOME \
# #         --libType A \
# #         --alignments $DATA/$SAMPLE/"$SAMPLE"_Aligned.toTranscriptome.out.bam \
# #         -p 10 \
# #         -g $ANNOTATION \
# #         -o $RESULT/$SAMPLE


# # Mapping-based

# declare -A SAMPLES=(
#     ['day0_1']='OJ63'
#     ['day0_2']='OJ64'
#     ['day0_3']='OJ65'
#     ['day0_4']='OJ66'
#     ['day0_5']='OJ67'
#     ['day1_1']='OJ68'
#     ['day1_2']='OJ69'
#     ['day1_3']='OJ70'
#     ['day1_4']='OJ71'
#     ['day1_5']='OJ72'
#     ['day2_1']='OJ73'
#     ['day2_2']='OJ74'
#     ['day2_3']='OJ75'
#     ['day2_4']='OJ76'
#     ['day2_5']='OJ77'
#     ['day3_1']='OJ79'
#     ['day3_2']='OJ80'
#     ['day3_3']='OJ81'
#     ['day3_4']='OJ82'
#     ['day3_5']='OJ83'
#     ['day5_1']='OJ85'
#     ['day5_2']='OJ86'
#     ['day5_4']='OJ87'
#     ['day5_5']='OJ88'
#     ['day5_6']='OJ89'
# )

# RAWSAMPLE=${SAMPLES["$SAMPLE"]}

# fa1=/project/owlmayerUnpublishedData/srRNAseq_differentiation_NGN3/rawData/"$RAWSAMPLE"/"$RAWSAMPLE"_R1.fastq.gz
# fa2=/project/owlmayerUnpublishedData/srRNAseq_differentiation_NGN3/rawData/"$RAWSAMPLE"/"$RAWSAMPLE"_R2.fastq.gz

# # echo $fa1
# salmon quant \
#         -i $INDEX \
#         --libType A \
#         -1 $fa1 \
#         -2 $fa2 \
#         -p 10 \
#         --validateMappings \
#         -o $RESULT/$SAMPLE