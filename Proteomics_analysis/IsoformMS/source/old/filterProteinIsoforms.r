list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))

#mode="all."
#mode="razor.unique."
mode="unique."

prefix="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/de-novo/LFQ."
MS.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/proteinGroups.csv"
annotation.info.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/de-novo/annotation.info.txt"

args <- commandArgs(trailingOnly = TRUE) 

MS.path = args[1]
annotation.info.path= args[2]
mode = args[3]
prefix= args[4]

ANNOTATION=read.table(annotation.info.path,header=T)
#TABLE_ALL=read.table(MS.path,sep="\t", header=T)
TABLE=read.table(MS.path,sep="\t", header=T)[,c("Protein.IDs","Peptide.counts..all.",
                                                paste("Peptide.counts..",mode,sep=""),
                                                "Number.of.proteins",
                                                "LFQ.intensity.Day0.B","LFQ.intensity.Day0.C","LFQ.intensity.Day0.D",
                                                "LFQ.intensity.Day1.B","LFQ.intensity.Day1.C","LFQ.intensity.Day1.D",
                                                "LFQ.intensity.Day2.B","LFQ.intensity.Day2.C","LFQ.intensity.Day2.D",
                                                "LFQ.intensity.Day3.B","LFQ.intensity.Day3.C","LFQ.intensity.Day3.D",              
                                                "LFQ.intensity.Day4.B","LFQ.intensity.Day4.C","LFQ.intensity.Day4.D",              
                                                "LFQ.intensity.Day5.B.C","LFQ.intensity.Day5.D.E","LFQ.intensity.Day5.F.G")]

TABLE=TABLE[rowSums(TABLE[,c("LFQ.intensity.Day0.B","LFQ.intensity.Day0.C","LFQ.intensity.Day0.D",
                             "LFQ.intensity.Day1.B","LFQ.intensity.Day1.C","LFQ.intensity.Day1.D",
                             "LFQ.intensity.Day2.B","LFQ.intensity.Day2.C","LFQ.intensity.Day2.D",
                             "LFQ.intensity.Day3.B","LFQ.intensity.Day3.C","LFQ.intensity.Day3.D",              
                             "LFQ.intensity.Day4.B","LFQ.intensity.Day4.C","LFQ.intensity.Day4.D",              
                             "LFQ.intensity.Day5.B.C","LFQ.intensity.Day5.D.E","LFQ.intensity.Day5.F.G")])!=0,]

selection=(TABLE[,"Peptide.counts..all."]==TABLE[,paste("Peptide.counts..",mode,sep="")]) |
  TABLE$Number.of.proteins==1
TABLE=as.data.table(TABLE[selection,])

TEMP_TABLE=TABLE[,.(transcript_id=unlist(strsplit(Protein.IDs,split=";"))),by="Protein.IDs"]
TEMP_TABLE=merge(TEMP_TABLE,ANNOTATION, by="transcript_id")
TEMP_TABLE=TEMP_TABLE[,.(transcript_id=transcript_id,
                         gene_id=paste(unique(gene_id),collapse=";"),
                         gene_name=gene_name
                         ),by="Protein.IDs"]

TEMP_TABLE=merge(TEMP_TABLE,TEMP_TABLE[,.(gene_name_unique=paste(unique(gene_name),collapse=";")),by="gene_id"],by="gene_id")
TEMP_TABLE=TEMP_TABLE[,.(Gene.IDs=paste(unique(gene_id),collapse=";"),
                         Gene.NAMEs=unique(gene_name_unique)),by="Protein.IDs"]

TABLE=merge(TABLE,TEMP_TABLE,by="Protein.IDs")

genes_of_interest=names(table(TABLE$Gene.IDs)[table(TABLE$Gene.IDs)>1])

TABLE=TABLE[TABLE$Gene.IDs%in%genes_of_interest,]
write.table(TABLE, col.names = T, row.names = F, sep = "\t", quote = F, file=paste(prefix,"tab",sep=""))
