#!/usr/bin/env bash

SOURCE=${BASH_SOURCE[0]}
DIR=$( dirname "$SOURCE" )

nanopore_gtf=Data/GRCh38_p12/nanopore.SQANTI3.gtf

out=Result/IsoformMS_SQANTI3_customProteome
mkdir -p $out

gtf_info=$out/annotation.info.txt

gtf=$nanopore_gtf

MS="Data/MS_NeuroDifferentation_NGN3/FlashLFQ/BayesianFoldChangeAnalysis.all.tsv"
LFQ="Data/MS_NeuroDifferentation_NGN3/FlashLFQ/FashLFQ.txt"
design="Data/MS_NeuroDifferentation_NGN3/ExperimentalDesign.continuous.txt"

cat <(cat Data/MS_NeuroDifferentation_NGN3/FlashLFQ/*.tabular | head -n 1) \
<(cat Data/MS_NeuroDifferentation_NGN3/FlashLFQ/*.tabular | grep -vP "^Protein Group\t") > $MS

ONT_RLE=Data/ONT-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt
Illumina_RLE=Data/RNA-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt

if [ ! -f $out/RNA-seq.txt ];then
	#This is the correct one, wait for ONT-seq re-quantification
	cat <(cat <(cat <(cat <(awk -v OFS="\t" -F"\t" '{ print $0" Comparison" }' Data/RNA-seq_NeuroDifferentation_NGN3/isoform/day1_background_day0_DEseq_all_res.csv \
	| head -n +1) \
	<(awk -v OFS="\t" -F"\t" '{ print $0" Day1_vs_Day0" }' Data/RNA-seq_NeuroDifferentation_NGN3/isoform/day1_background_day0_DEseq_all_res.csv | tail -n +2)) \
	<(awk -v OFS="\t" -F"\t" '{ print $0" Day2_vs_Day1" }' Data/RNA-seq_NeuroDifferentation_NGN3/isoform/day2_background_day1_DEseq_all_res.csv | tail -n +2)) \
	<(awk -v OFS="\t" -F"\t" '{ print $0" Day3_vs_Day2" }' Data/RNA-seq_NeuroDifferentation_NGN3/isoform/day3_background_day2_DEseq_all_res.csv | tail -n +2)) \
	<(awk -v OFS="\t" -F"\t" '{ print $0" Day5_vs_Day3" }' Data/RNA-seq_NeuroDifferentation_NGN3/isoform/day5_background_day3_DEseq_all_res.csv | tail -n +2) \
	> $out/RNA-seq.txt
fi

if [ ! -f $out/ONT-seq.txt ];then
	#This is the correct one, wait for ONT-seq re-quantification
	cat <(cat <(awk -v OFS="\t" -F"\t" '{ print $0" Comparison" }' Data/ONT-seq_NeuroDifferentation_NGN3/isoform/day3_vs_day0_DEseq_all_res.csv \
	| head -n +1) \
	<(awk -v OFS="\t" -F"\t" '{ print $0" Day3_vs_Day0" }' Data/ONT-seq_NeuroDifferentation_NGN3/isoform/day3_vs_day0_DEseq_all_res.csv | tail -n +2)) \
	<(awk -v OFS="\t" -F"\t" '{ print $0" Day5_vs_Day3" }' Data/ONT-seq_NeuroDifferentation_NGN3/isoform/day5_vs_day3_DEseq_all_res.csv | tail -n +2) \
	> $out/ONT-seq.txt
fi

if [ ! -f $gtf_info ];then
	cat <(echo "gene_id	transcript_id	gene_name") \
	<(	
	paste <(paste <(grep -P "\ttranscript\t" $gtf | cut -f9 | grep -Po "gene_id .*"| cut -f1 -d';' | cut -f2 -d'"') \
	<(grep -P "\ttranscript\t" $gtf | cut -f9 | grep -Po "transcript_id .*"| cut -f1 -d';' | cut -f2 -d'"')) \
	<(grep -P "\ttranscript\t" $gtf | cut -f9 |cut -f3 -d';' | cut -f2 -d'"') \
	) \
	> $gtf_info
fi

ratio=$out/MSratio.atLeastTwoProteinIsoforms.tab
fdr=$out/MSsd.atLeastTwoProteinIsoforms.tab
if [ ! -f $ratio ] || [ ! -f $fdr ];then
		Rscript $DIR/source/prepareFlashLFQoutput.r $MS $gtf_info 14 atLeastTwoProteinIsoforms $ratio $fdr \
		"Day1_vs_Day0;Day2_vs_Day1;Day3_vs_Day2;Day4_vs_Day3;Day5_vs_Day4"
fi

ratio_illumina=$out/MSratio.atLeastTwoProteinIsoforms.Illumina.tab
fdr_illumina=$out/MSsd.atLeastTwoProteinIsoforms.Illumina.tab

if [ ! -f $ratio_illumina ] || [ ! -f $fdr_illumina ];then
		Rscript $DIR/source/prepareFlashLFQoutput.r $MS $gtf_info 14 atLeastTwoProteinIsoforms $ratio_illumina $fdr_illumina \
		"Day1_vs_Day0;Day2_vs_Day1;Day3_vs_Day2;Day5_vs_Day3"
fi

ratio_ont=$out/MSratio.atLeastTwoProteinIsoforms.ont.tab
fdr_ont=$out/MSsd.atLeastTwoProteinIsoforms.ont.tab

if [ ! -f $ratio_ont ] || [ ! -f $fdr_ont ];then
		Rscript $DIR/source/prepareFlashLFQoutput.r $MS $gtf_info 14 atLeastTwoProteinIsoforms $ratio_ont $fdr_ont \
		"Day3_vs_Day0;Day5_vs_Day3"
fi

# get number of overlap transcripts for each transcript
# idea: non-overlapping transcripts should be removed
if [ ! -f $out/overlap.tab ];then
	bedtools=`which bedtools`

	transcriptIDs=`cut -f1 $ratio | tail -n +2 | sed "s/;/\n/g" | sort | uniq | \
	awk -v OFS="\t" -F"\t" '{ print "transcript_id \""$0"\"" }' | tr '\n' '|'`

	cut -f1 $ratio | tail -n +2 | sed "s/;/\n/g" | sort | uniq | wc -l

	grep -P "\ttranscript\t" $gtf | grep -P "${transcriptIDs::-1}" > $out/temp.gtf 
	$bedtools intersect -c -s -a $out/temp.gtf -b $out/temp.gtf > $out/overlap.gtf

	paste <(grep -P "transcript_id \"*\"" $out/overlap.gtf | cut -f2 -d'"') \
	<(cut -f10 $out/overlap.gtf) > $out/overlap.tab
	
	rm -f $out/temp.gtf 
	rm -f $out/overlap.gtf
fi

if [ ! -f $out/trendChanges.heatmap.pdf ];then
	Rscript $DIR/source/plotTrendDifference.r $ratio $LFQ $out/overlap.tab 5 1 $design $out/trendChanges.pdf $out/trendChanges.heatmap.pdf
fi

if [ ! -f $out/correlation.Illumina.pdf ];then
	Rscript $DIR/source/plotCorrelation.r $ratio_illumina $out/overlap.tab $out/RNA-seq.txt $out/correlation.Illumina.pdf
fi

if [ ! -f $out/correlation.ONT.pdf ];then
	Rscript $DIR/source/plotCorrelation.r $ratio_ont $out/overlap.tab $out/ONT-seq.txt $out/correlation.ONT.pdf
fi

if [ ! -f $out/integrativePlot.pdf ];then
	Rscript $DIR/source/plotIntegrativeVisualization.r $ratio $LFQ $out/overlap.tab $design $Illumina_RLE $ONT_RLE $out/integrativePlot.pdf
fi

