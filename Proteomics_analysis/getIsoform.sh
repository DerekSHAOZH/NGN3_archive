#!/usr/bin/env bash

nanopore_gtf=../Data/GRCh38_p12/nanopore.SQANTI3.gtf

outDir=../Result/IsoformMS_SQANTI3
mkdir -p $outDir


# get de-novo translation information	
if [ ! -f $outDir/translation.de-novo.tab ];then
	#get coding sequence fron all de-novo transcripts
	if [ ! -f $outDir/nanopore.transcript.fa ];then
		gffread=`which gffread`
		$gffread \
		-w $outDir/nanopore.transcript.fa \
		-g ../Data/GRCh38_p12/GRCh38.p12.primary_assembly.genome.fa \
		$nanopore_gtf
	fi
	
	#get AA sequence from all de-novo transcripts 
	if [ ! -f $outDir/nanopore.translation.fa ];then
		python source/translation_to_protein.py $outDir/nanopore.transcript.fa $outDir/nanopore.translation.fa
	fi
	
	if [ ! -f $outDir/translation.de-novo.coord ];then
		awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' \
		$outDir/nanopore.translation.fa | sed -r "s/^>//g" | \
		sort -k1,1 > $outDir/translation.de-novo.coord
	fi
	
	if [ ! -f $outDir/translation.de-novo.tab ];then
		cut -f1,4 $outDir/translation.de-novo.coord > $outDir/translation.de-novo.tab
	fi
	
	rm -f $outDir/nanopore.transcript.fa
	rm -f $outDir/nanopore.translation.fa	
fi

mode="de-novo" 

out=$outDir/$mode
mkdir -p $out
	
gtf_info=$out/annotation.info.txt

gtf=$nanopore_gtf

	
if [ ! -f $gtf_info ];then
	cat <(echo "gene_id	transcript_id") \
	<(	
	paste <(grep -P "\ttranscript\t" $gtf | cut -f9 | grep -Po "gene_id .*"| cut -f1 -d';' | cut -f2 -d'"') \
	<(grep -P "\ttranscript\t" $gtf | cut -f9 | grep -Po "transcript_id .*"| cut -f1 -d';' | cut -f2 -d'"')
	) \
	> $gtf_info
fi

if [ ! -f $out/PeptideToTranslationMapping.txt ];then
	# run on computing cluster with at least 70 nodes such as holidayincambodia or elcattivo
	Rscript source/searchPeptides.r \
	../Data/MS_NeuroDifferentation_NGN3/evidence.txt \
	$gtf_info \
	$outDir/translation.$mode.tab \
	$out/PeptideToTranslationMapping.txt > $out/searchPeptides.log
fi
	
if [ ! -f $out/PeptidesMappingDistinctTranscriptGroups.txt ];then
	Rscript source/getPeptidesMappingDistinctTranscriptGroups.r \
	../Data/MS_NeuroDifferentation_NGN3/evidence.txt \
	$out/PeptideToTranslationMapping.txt \
	$out/PeptidesMappingDistinctTranscriptGroups.txt
fi
	
# get number of overlap transcripts for each transcript
# idea: non-overlapping transcripts should be removed
if [ ! -f $out/transcriptsWithUniquePeptide.overlap.tab ];then
	bedtools=`which bedtools`
	
	transcriptIDs=`cut -f10 $out/PeptidesMappingDistinctTranscriptGroups.txt | tail -n +2 | sort | uniq | \
	sed "s/;/\n/g" | sed "s/|/\n/g" | awk -v OFS="\t" -F"\t" '{ print "transcript_id \""$0"\"" }' | tr '\n' '|'`

	grep -P "${transcriptIDs::-1}" $gtf | grep -P "\ttranscript\t" > $out/transcriptsWithUniquePeptide.gtf 
	$bedtools intersect -c -s -a $out/transcriptsWithUniquePeptide.gtf -b $out/transcriptsWithUniquePeptide.gtf \
	> $out/transcriptsWithUniquePeptide.overlap.gtf

	paste <(grep -P "transcript_id \"*\"" $out/transcriptsWithUniquePeptide.overlap.gtf | cut -f2 -d'"') \
	<(cut -f10 $out/transcriptsWithUniquePeptide.overlap.gtf) > $out/transcriptsWithUniquePeptide.overlap.tab
	
	rm -f $out/transcriptsWithUniquePeptide.overlap.gtf
	rm -f $out/transcriptsWithUniquePeptide.gtf
fi
	
if [ ! -f $out/IsoformIntensities.txt ] || [ ! -f $out/pdf/Isoform.Intensities.norm1.pdf ];then
	mkdir -p $out/pdf
	Rscript source/plotIsoformChanges.r \
	$out/PeptidesMappingDistinctTranscriptGroups.txt \
	$gtf \
	$out/transcriptsWithUniquePeptide.overlap.tab \
	2 2 2 \
	$out/IsoformIntensities.txt \
	$out/pdf/Isoform.Intensities.norm1.pdf > $out/pdf/plotIsoformTrendChanges.log
fi
	
if [ ! -f $out/IsoformTrendchanges.txt ] || [ ! -f $out/pdf/Isoform.Trendchanges.logFC.pdf ];then
	Rscript source/plotIsoformTrendChange.r \
	$out/IsoformIntensities.txt \
	$out/IsoformTrendchanges.txt \
	$out/pdf/Isoform.Trendchanges.logFC.pdf \
	$out/pdf/Isoform.Trendchanges.logFC.mean.pdf
fi
	
if [ ! -f $out/pdf/Isoform.Trendchanges.heatmap.noChangesRemove.pdf ];then
	mkdir -p $out/ids
	Rscript source/plotIsoformTrendChangeHeatmap.r $out/IsoformTrendchanges.txt \
	$out/pdf/Isoform.Trendchanges.heatmap.pdf $out/pdf/Isoform.Trendchanges.heatmap.noChangesRemove.pdf \
	$out/ids/
fi

	

