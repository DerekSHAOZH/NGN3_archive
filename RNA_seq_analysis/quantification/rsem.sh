##!/usr/bin/env python
#! /usr/bin/bash

GENOME=/project/Neurodifferentiation_System/GeneralResources/GRCh38.p12.primary_assembly.genome.fa
WORKING=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Scripts
DATA=/project/owlmayerTemporary/Derek/NGN3_project/RNA_alignGenome/alignment
rsem=/project/owlmayer/Applications/RSEM-master

# File=nanopore_corrected
File=nanopore.filtered


## rsem-prepare-reference




# # TRANSCRIPTOME=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/"$File".gtf
# TRANSCRIPTOME=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/"$File".gtf


# RESULT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Data/rsem_"$File"_index


# mkdir -p $WORKING/mxq_log/AlignToOntTrans
# mkdir -p $RESULT
 
# mxqsub -w $WORKING -o $WORKING/mxq_log/AlignToOntTrans/rsem_"$File"_index.log -j 10 -m 50G -t 10h \
# 			$rsem/rsem-prepare-reference --gtf $TRANSCRIPTOME \
#                --star \
#                $GENOME \
#                $RESULT/$File


## rsem-calculate-expression
reference=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Data/rsem_"$File"_index/"$File"

declare -A SAMPLES=(
    ['day0_1']='OJ63'
    ['day0_2']='OJ64'
    ['day0_3']='OJ65'
    ['day0_4']='OJ66'
    ['day0_5']='OJ67'
    ['day1_1']='OJ68'
    ['day1_2']='OJ69'
    ['day1_3']='OJ70'
    ['day1_4']='OJ71'
    ['day1_5']='OJ72'
    ['day2_1']='OJ73'
    ['day2_2']='OJ74'
    ['day2_3']='OJ75'
    ['day2_4']='OJ76'
    ['day2_5']='OJ77'
    ['day3_1']='OJ79'
    ['day3_2']='OJ80'
    ['day3_3']='OJ81'
    ['day3_4']='OJ82'
    ['day3_5']='OJ83'
    ['day5_1']='OJ85'
    ['day5_2']='OJ86'
    ['day5_4']='OJ87'
    ['day5_5']='OJ88'
    ['day5_6']='OJ89'
)



for subdir in $DATA/day*/
do
    # filename=$(basename -- "$f")
    # extension="${filename##*.}"
    # filename="${filename%.*}"

    subdir_name=$(basename "$subdir")
    # echo "Subdirectory name: $subdir_name"

    pathToOutputFiles=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/RSEM/"$subdir_name"

    mkdir -p $pathToOutputFiles

    RAWSAMPLE=${SAMPLES["$subdir_name"]}

    fa1=/project/owlmayerUnpublishedData/srRNAseq_differentiation_NGN3/rawData/"$RAWSAMPLE"/"$RAWSAMPLE"_R1.fastq.gz
    fa2=/project/owlmayerUnpublishedData/srRNAseq_differentiation_NGN3/rawData/"$RAWSAMPLE"/"$RAWSAMPLE"_R2.fastq.gz
    

    # mxqsub -w $pathToOutputFiles -o $pathToOutputFiles/$subdir_name.log -j 10 -m 100G -t 15h \
    # $rsem/rsem-calculate-expression -p 10 --paired-end --star --star-gzipped-read-file --star-output-genome-bam $fa1 $fa2 $reference $subdir_name #\
    # samtools index $subdir_name.STAR.genome.bam $subdir_name.STAR.genome.bam.bai


    # zcat $fa1 > $pathToOutputFiles/p.fastq \
    # zcat $fa2 > $pathToOutputFiles/n.fastq \
    # $rsem/rsem-calculate-expression -p 10 --paired-end $pathToOutputFiles/p.fastq $pathToOutputFiles/n.fastq $reference $subdir_name \
    rm $pathToOutputFiles/*.fastq 

    # $rsem/rsem-calculate-expression -p 10 --paired-end <($fa1) <($fa2) $reference $subdir_name \

    
   

    # echo $reference
    # echo "$subdir_name"
done





