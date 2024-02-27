#!/usr/bin/env bash

SOURCE=${BASH_SOURCE[0]}
DIR=$( dirname "$SOURCE" )

outDir=Result/RBPs

mkdir -p $outDir

if [ ! -f $outDir/significant.txt ];then
	cat <(cat Data/ONT-seq_NeuroDifferentation_NGN3/gene/*_vs_*_DEseq_all_res.csv | awk -v OFS="," -F"," '{ if ( $6 <= 0.05) print $7}' | sort | uniq) \
	<(cat Data/RNA-seq_NeuroDifferentation_NGN3/gene/day[3,5]_background_day[0,3]_DEseq_all_res.csv | \
	awk -v OFS=" " -F" " '{ if ( $6 <= 0.05) print $7}' | sort | 	uniq) | \
	sort | uniq > $outDir/significant.txt
fi

if [ ! -f $outDir/RBP_splicing.ONT.pdf ] || [ ! -f $outDir/RBP_splicing.RNA.pdf ];then
	Rscript $DIR/source/plotGOHeatmap.integrative.r "GO:0008380;GO:0003723" \
	Data/ONT-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt \
	Data/RNA-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt \
	$outDir/significant.txt \
	4 $outDir/RBP_splicing.txt $outDir/RBP_splicing.ONT.pdf $outDir/RBP_splicing.RNA.pdf \
	$outDir/cluster.txt $outDir/reference.cluster.txt
fi

grep -P "\t1$" $outDir/cluster.txt | cut -f1 > $outDir/cluster.1.txt
grep -P "\t2$" $outDir/cluster.txt | cut -f1 > $outDir/cluster.2.txt
grep -P "\t3$" $outDir/cluster.txt | cut -f1 > $outDir/cluster.3.txt
grep -P "\t4$" $outDir/cluster.txt | cut -f1 > $outDir/cluster.4.txt

if [ -f Data/GOfiles/cluster_1.GOslim_Biological_Process.json ] && [ ! -f $outDir/cluster_1.GOslim_Biological_Process.pdf ];then
	Rscript $DIR/source/plotGO.json.r Data/GOfiles/cluster_1.GOslim_Biological_Process.json 0 10 "+" $outDir/cluster_1.GOslim_Biological_Process.pdf
fi

#if [ -f Data/GOfiles/cluster_3.GOslim_Biological_Process.json ] && [ ! -f $outDir/cluster_3.GOslim_Biological_Process.pdf ];then
	Rscript $DIR/source/plotGO.json.r Data/GOfiles/cluster_3.GOslim_Biological_Process.json 0 10 "+" $outDir/cluster_3.GOslim_Biological_Process.pdf
#fi

