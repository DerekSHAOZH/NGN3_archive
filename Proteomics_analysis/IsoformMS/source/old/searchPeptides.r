list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


suppressPackageStartupMessages(library(data.table))


MS.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/evidence.txt"
annotation.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS/de-novo/annotation.info.txt"
AS.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS/translation.de-novo.tab"
out="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS/de-novo/PeptideToTranslationMapping.txt"

args <- commandArgs(trailingOnly = TRUE) 

MS.path = args[1]
annotation.path = args[2]
AS.path = args[3]
out= args[4]

#prepare MS table
MS=read.table(MS.path, sep="\t", header=T)

MS=MS[,c("Sequence", "Modifications","Intensity", "Gene.names", "Experiment", "Peptide.ID")]
MS=MS[MS$Modifications=="Unmodified",]
MS$Day = do.call(rbind,strsplit(MS$Experiment,split="-"))[,1]
MS=MS[!is.na(MS$Intensity),]

#prepare amino acid translations
AS=fread(AS.path, sep="\t", header=F, col.names=c("transcript_id","sequence"))

# prepare annotation reference
annotation = fread(annotation.path,sep="\t", header=T)

#63,738 transcripts
TABLE = merge(annotation, AS, by="transcript_id")
print("# reference transcripts")
print(nrow(TABLE))

TABLE = TABLE[TABLE$sequence!="M",]
print("# putative translated transcripts")
print(nrow(TABLE))

# 56,374 unique protein sequences
TABLE=TABLE[,.(transcript_id=paste(transcript_id,collapse="|"),
              gene_id=paste(unique(gene_id),collapse="|")
              ), by="sequence"]
print("# unique AA sequences")
print(nrow(TABLE))

#identify for each peptide occurence in reference translation
Peptides=unique(MS[,c("Sequence","Peptide.ID")])
print("# Peptides")
print(nrow(Peptides))

get_transcriptID = function(X,Y){
  return(paste(sort(Y$transcript_id[grepl(X[1],Y$sequence, fixed=T)]), collapse=";"))
}

get_geneID= function(X,Y){
  return(paste(sort(unique(Y$gene_id[grepl(X[1],Y$sequence, fixed=T)])), collapse=";"))
}

#this block takes a few minutes, better run on interactive servers, e-g- holidayincambodia or elcattivo
library(future.apply)
plan(multisession, workers = 70) 
Peptides$transcript_id=future_apply(Peptides,1,get_transcriptID,Y=TABLE)
Peptides$gene_id=future_apply(Peptides,1,get_geneID,Y=TABLE)

print("# assigned peptides")
print(sum(Peptides$transcript_id!=""))

write.table(Peptides, col.names = T, row.names = F, sep = "\t", quote = F, file = out)
