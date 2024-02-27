#!/usr/bin/env bash

SOURCE=${BASH_SOURCE[0]}
DIR=$( dirname "$SOURCE" )

# get relevant data
mkdir -p Data/GRCh38_p12
mkdir -p Data/MS_NeuroDifferentation_NGN3
mkdir -p Data/ONT-seq_NeuroDifferentation_NGN3/isoform
mkdir -p Data/ONT-seq_NeuroDifferentation_NGN3/gene
mkdir -p Data/RNA-seq_NeuroDifferentation_NGN3/isoform
mkdir -p Data/RNA-seq_NeuroDifferentation_NGN3/gene

############# move this to data preparation

ln -sf /project/Neurodifferentiation_System/GeneralResources/gencode.v32.primary_assembly.annotation.gtf Data/GRCh38_p12/gencode.v32.primary_assembly.annotation.gtf

if [ ! -f Data/GRCh38_p12/nanopore.SQANTI3.proteome.fasta ];then
	awk -v OFS="\t" -F"\t" '{ if ( $0 ~ />/ ) print $1; else print $0}' \
	/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore.filtered.faa > Data/GRCh38_p12/nanopore.SQANTI3.proteome.fasta
fi


##if [ ! -f Data/GRCh38_p12/nanopore.SQANTI3.proteome.lookup.txt ];then
	join -t $'\t' -j 1 <(awk -v OFS="\t" -F"\t" '{ if ( $0 ~ />/ ) print substr($1,2);}' Data/GRCh38_p12/nanopore.SQANTI3.proteome.fasta | sort) \
	<(cut -f1,33,34,47 /project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_RulesFilter_result_classification.txt | \
	sort -k1,1) > Data/GRCh38_p12/nanopore.SQANTI3.proteome.lookup.txt
#fi

# SQUANTI filtering does not contain gene names, take gffcompare output and select only transcripts from SQUANTI filtering
if [ ! -f Data/GRCh38_p12/nanopore.SQANTI3.gtf ];then
	cp /project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore.filtered.gtf Data/GRCh38_p12/nanopore.SQANTI3.temp.gtf
	
	Rscript $DIR/source/filterGTF.r Data/GRCh38_p12/nanopore.SQANTI3.temp.gtf /project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/GffCompare/nanopore.combined.filt.gtf Data/GRCh38_p12/nanopore.SQANTI3.temp2.gtf
	
	grep -Pv "^#" Data/GRCh38_p12/nanopore.SQANTI3.temp2.gtf > Data/GRCh38_p12/nanopore.SQANTI3.gtf
	rm Data/GRCh38_p12/nanopore.SQANTI3.temp*
fi

if [ ! -f Data/GRCh38_p12/GRCh38.p12.primary_assembly.genome.fa ];then
	ln -fs /project/owlmayerTemporary/Assembly/Human/GRCh38_p12/GRCh38.p12.primary_assembly.genome.fa Data/GRCh38_p12/GRCh38.p12.primary_assembly.genome.fa
fi

if [ ! -f Data/ONT-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt ];then
	sed "s/,/ /g" /project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Results/DTE_analysis/count_transcript_deseq2norm.csv | \
	sed "s/day/Day/g" > \
	Data/ONT-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt
fi

if [ ! -f Data/RNA-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt ];then
	sed "s/,/ /g" /project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/DTE_analysis_with_RSEM/count_transcript_deseq2norm.csv | sed "s/day/Day/g" > \
	Data/RNA-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt
fi
	
if [ ! -f Data/ONT-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt ];then
	sed "s/,/ /g" /project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/DGE_analysis/count_deseq2norm.csv | \
	sed "s/day/ONT_day/g" > \
	Data/ONT-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt
fi

if [ ! -f Data/RNA-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt ];then
	sed "s/day/RNA_day/g" /project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/counts_normalized.txt > \
	Data/RNA-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt
fi

for file in `ls /project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/*/DEseq_all_res.csv`
do
	name=`basename $(dirname $file)`
	if [ ! -f Data/RNA-seq_NeuroDifferentation_NGN3/gene/${name}_DEseq_all_res.csv ];then
		cp $file Data/RNA-seq_NeuroDifferentation_NGN3/gene/${name}_DEseq_all_res.csv
	fi
done

for file in `ls /project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/DTE_analysis_with_RSEM/*/DEseq_all_res.csv`
do
	name=`basename $(dirname $file)`
	if [ ! -f Data/RNA-seq_NeuroDifferentation_NGN3/isoform/${name}_DEseq_all_res.csv ];then
		cp $file Data/RNA-seq_NeuroDifferentation_NGN3/isoform/${name}_DEseq_all_res.csv
	fi
done

for file in `ls /project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Results/DTE_analysis/*/DEseq_all_res.csv`
do
	name=`basename $(dirname $file)`
	if [ ! -f Data/ONT-seq_NeuroDifferentation_NGN3/isoform/${name}_DEseq_all_res.csv ];then
		sed "s/,/ /g" $file > Data/ONT-seq_NeuroDifferentation_NGN3/isoform/${name}_DEseq_all_res.csv
	fi
done

for file in `ls /project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/DGE_analysis/*/DEseq_all_res.csv`
do
	name=`basename $(dirname $file)`
	if [ ! -f Data/ONT-seq_NeuroDifferentation_NGN3/gene/${name}_DEseq_all_res.csv ];then
		cp $file Data/ONT-seq_NeuroDifferentation_NGN3/gene/${name}_DEseq_all_res.csv
	fi
done

