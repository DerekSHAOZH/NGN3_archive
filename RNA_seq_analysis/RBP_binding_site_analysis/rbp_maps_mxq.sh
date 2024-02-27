##!/usr/bin/env python
#! /usr/bin/bash

# source /home/shao/Tools/miniconda3/etc/profile.d/conda.sh
# conda activate rbp-maps





###TO-CHANGE
# EXPERIMENT=$1
# COMPARISON=$2

WORKING=/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_RBP
DATA=/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Data/CLIPdb_RBP
SCRIPT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Scripts/rbp_maps_run.sh
comparison_list=("day1_background_day0" "day2_background_day1" "day3_background_day2" "day5_background_day3")


# for f in $WORKING/*.bb
for f in $DATA/*.bb
do
    filename=$(basename -- "$f")
    extension="${filename##*.}"
    filename="${filename%.*}"
    
    for comparison in "${comparison_list[@]}"; do
    	OUTPUT=$WORKING/"$filename"/"$comparison"

    	# echo $OUTPUT
    	mkdir -p $OUTPUT

    	mxqsub -w $OUTPUT -o $OUTPUT/rbp_maps.log -j 10 -m 50G -t 10h sh $SCRIPT $filename $comparison

	done



done


# RMATS=/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitches_RBPs/rMAPS/Results/"$COMPARISON"/SE

# mkdir -p $OUTPUT

# PEAK=$WORKING/peaks_final.bed

# UP_RMATS=$RMATS/up.SE.MATS.JC.txt
# DN_RMATS=$RMATS/dn.SE.MATS.JC.txt
# BG_RMATS=$RMATS/bg.SE.MATS.JC.txt

# # RMATS_NR=$DATA/SE.MATS.JC.nr.txt
# # mkdir -p $OUTPUT
# # mkdir -p $LOG


# # subset_jxc -i $RMATS \
# # -o $DATA/SE.MATS.JC.nr.txt \
# # -e se


# plot_map --peak $PEAK --annotations $UP_RMATS $DN_RMATS $BG_RMATS \
# --annotation_type rmats rmats rmats --output $DATA/test.png \
# --event se \
# --normalization_level 0 --testnums 0 1 --bgnum 2 --sigtest fisher \
# --intron_offset 250

# for f in $DATA/*.Wig
# do
#     filename=$(basename -- "$f")
#     extension="${filename##*.}"
#     filename="${filename%.*}"

#     # bam=$OUTPUT/"$filename".bam 		
#     # sbam=$OUTPUT/"$filename".sorted.bam
#     # sam=$OUTPUT/"$filename".sam

#     # mxqsub -w $WORKING -o $LOG/"$filename"_map2transcriptome.log -j 10 -m 100G -t 10h \
#     # rm $DATA/"$filename".bed
#     wig2bed --zero-indexed < $f | bedmap --max $PEAK - > $DATA/"$filename"_height.bed;

#     # mxqsub -w $WORKING -j 10 -m 100G -t 10h \
#     # minimap2 -ax map-ont -t 10 -p $psec -N $msec $INDEX $f | samtools view -Sb > $bam
#         # samtools sort -@ 10 $bam -o $sbam;
#         # samtools index $sbam;



#     # subdir_name=$(basename "$subdir")
#     # # echo "Subdirectory name: $subdir_name"

#     # pathToOutputFiles=/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/RSEM/"$subdir_name"

#     # mkdir -p $pathToOutputFiles

#     # RAWSAMPLE=${SAMPLES["$subdir_name"]}

#     # fa1=/project/owlmayerUnpublishedData/srRNAseq_differentiation_NGN3/rawData/"$RAWSAMPLE"/"$RAWSAMPLE"_R1.fastq.gz
#     # fa2=/project/owlmayerUnpublishedData/srRNAseq_differentiation_NGN3/rawData/"$RAWSAMPLE"/"$RAWSAMPLE"_R2.fastq.gz
    

#     # # mxqsub -w $pathToOutputFiles -o $pathToOutputFiles/$subdir_name.log -j 10 -m 100G -t 15h \
#     # # $rsem/rsem-calculate-expression -p 10 --paired-end --star --star-gzipped-read-file --star-output-genome-bam $fa1 $fa2 $reference $subdir_name #\
#     # # samtools index $subdir_name.STAR.genome.bam $subdir_name.STAR.genome.bam.bai


#     # # zcat $fa1 > $pathToOutputFiles/p.fastq \
#     # # zcat $fa2 > $pathToOutputFiles/n.fastq \
#     # # $rsem/rsem-calculate-expression -p 10 --paired-end $pathToOutputFiles/p.fastq $pathToOutputFiles/n.fastq $reference $subdir_name \
#     # rm $pathToOutputFiles/*.fastq 

#     # # $rsem/rsem-calculate-expression -p 10 --paired-end <($fa1) <($fa2) $reference $subdir_name \

    
   

#     # echo $reference
#     # echo $f
# done

# paste $PEAK $DATA/*height.bed > $DATA/answer.bed

   
