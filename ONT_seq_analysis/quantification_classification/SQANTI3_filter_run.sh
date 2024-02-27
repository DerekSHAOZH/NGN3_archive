##!/usr/bin/env python
#! /usr/bin/bash

source /home/shao/Tools/miniconda3/etc/profile.d/conda.sh
conda activate SQANTI3.env

export PYTHONPATH=$PYTHONPATH:~/Tools/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:~/Tools/cDNA_Cupcake/

echo "$(date): SQANTI3 filter starts"


TOOL=/home/shao/Tools/SQANTI3-5.1.1/sqanti3_filter.py
GTF=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_corrected.gtf
# REFGTF=/project/Neurodifferentiation_System/GeneralResources/gencode.v32.primary_assembly.annotation.gtf
# REFFASTA=/project/Neurodifferentiation_System/GeneralResources/GRCh38.p12.primary_assembly.genome.fa
FAA=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_corrected.faa
# CAGE=/home/shao/Tools/SQANTI3-5.1.1/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed
# POLYA=/home/shao/Tools/SQANTI3-5.1.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt
# POLYASITE=/home/shao/Tools/SQANTI3-5.1.1/data/polyA_site/atlas.clusters.2.0.GRCh38.96.bed
# COUNT=/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_sqanti3.txt
# SJ=/project/owlmayerTemporary/Derek/RNA_alignGenome/sr_Junction
RULE=/home/shao/Tools/SQANTI3-5.1.1/utilities/filter/filter_nanopore.json
CLASS=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_classification.txt
OUTPUT=/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter

python $TOOL rules --gtf $GTF --faa $FAA -o nanopore -d $OUTPUT -j $RULE $CLASS    




   
