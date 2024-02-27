#!/usr/bin/env bash


#########
TRANSLATION=/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/GRCh38_p12/nanopore.SQANTI3.proteome.lookup.txt
#pep=/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_noSingleExon/de-novo/PeptidesMappingDistinctTranscriptGroups.txt

gtf=/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/GRCh38_p12/nanopore.SQANTI3.gtf
peptide_col="#0d00ff"
background="#b3b3b3"
#background="remove"

outDir=/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/IsoformMS/processed/test
mkdir -p $outDir
######################
gene="PFN2"
gene_id="XLOC_032424"
mark=\
"TCONS_00083909/col1
TCONS_00083910/col2"

nT=`echo $mark | awk '{print NF}'`

if [ "$nT" == "2" ];then
	col1="#f8766d"
	col2="#00bfc4"
	markT=`echo $mark | sed "s/col1/$col1/g" | sed "s/col2/$col2/g"`
fi

if [ "$nT" == "3" ];then
	col1="#f8766d"
	col2="#00ba38"
	col3="#619cff"
	markT=`echo $mark | sed "s/col1/$col1/g" | sed "s/col2/$col2/g" | sed "s/col3/$col3/g"`
fi

if [ "$nT" == "4" ];then
	col1="#f8766d"
	col2="#7cae00"
	col3="#00bfc4"
	col4="#c77cff"
	markT=`echo $mark | sed "s/col1/$col1/g" | sed "s/col2/$col2/g" | sed "s/col3/$col3/g" | sed "s/col4/$col4/g"`
fi
echo $markT
############################

mkdir -p $outDir

grep -P "$gene_id" $gtf | gffread - -L -F -o $outDir/$gene.gff3
### ADD CDS

Rscript source/getCDS.r $outDir/$gene.gff3 $outDir/$gene.CDS.gff3 $TRANSLATION
cp $outDir/$gene.CDS.gff3 $outDir/$gene.CDS.temp.gff3


for row in $markT
do
	for transcript in `echo $row | cut -f1 -d'/' | sed "s/|/ /g" | sed "s/;/ /g"`
	do
		#transcript=`echo $row | cut -f1 -d'/'`
		col=`echo $row | cut -f2 -d'/'`
	
		echo $transcript
		echo $col
		awk -v OFS="\t" -F"\t" '{ if ( $0 ~ /'$transcript'/ ) print $0";color='$col'"; else print $0}' $outDir/$gene.CDS.temp.gff3 > $outDir/$gene.CDS.temp2.gff3 
		mv $outDir/$gene.CDS.temp2.gff3 $outDir/$gene.CDS.temp.gff3
	done
done

if [ "$background" == "remove" ];then
	awk -v OFS="\t" -F"\t" '{ if ( $0 ~ /color=/ ) print $0}' $outDir/$gene.CDS.temp.gff3 \
	> $outDir/$gene.CDS.mark.gff3
else
	awk -v OFS="\t" -F"\t" '{ if ( $0 ~ /color=/ ) print $0; else print $0";color='$background'"}' $outDir/$gene.CDS.temp.gff3 \
	> $outDir/$gene.CDS.mark.gff3
fi

rm -f $outDir/$gene.gff3
rm -f $outDir/$gene.CDS.gff3
rm -f $outDir/$gene.CDS.temp.gff3

#Rscript source/addUniquePeptide.r $outDir/$gene.CDS.mark.gff3 $outDir/$gene.CDS.mark.pep.gff3 $pep $de_novo
#rm -f $outDir/$gene.CDS.mark.gff3
