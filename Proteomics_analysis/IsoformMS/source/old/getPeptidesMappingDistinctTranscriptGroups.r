list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


suppressPackageStartupMessages(library(data.table))


peptides.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_noSingleExon/de-novo/PeptideToTranslationMapping.txt"
MS.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/evidence.txt"

args <- commandArgs(trailingOnly = TRUE) 

MS.path = args[1]
peptides.path = args[2]
out=args[3]

#prepare MS table
MS=read.table(MS.path, sep="\t", header=T)

MS=MS[,c("Sequence", "Modifications","Intensity", "Gene.names", "Experiment", "Peptide.ID")]
MS=MS[MS$Modifications=="Unmodified",]
MS$Day = do.call(rbind,strsplit(MS$Experiment,split="-"))[,1]
MS=MS[!is.na(MS$Intensity),]

#load peptides
Peptides = fread(peptides.path)

#remove unassigned peptides
Peptides = Peptides[Peptides$transcript_id!="",]

is.in = function(X,Y){
  return(Y%in%X)
}
  
help_internal = function(Y, is_unique){
  return(sum(Y %in% is_unique)==length(Y))
}

help = function(X){
  split= strsplit(X, split=";")
  lsplit=lengths(split)
  
  is_unique = names(table(unlist(split)))[table(unlist(split))==1]
  
  
  result=unlist(lapply(split, help_internal, is_unique=is_unique))
  result[lsplit==1]=T
  return(result)
}

#collect peptide IDs that fit same gene and transcript id
Peptides_cluster = Peptides[,.(Peptide.IDs=paste(unique(Peptide.ID),collapse=";")), by=c("gene_id","transcript_id")]

#identify peptide groups that map uniquely to non-overlapping transcript group
##### get if transcript group is smallest non overlapping group
Peptides_cluster = Peptides_cluster[,.(Peptide.IDs=Peptide.IDs, transcript_id = transcript_id, smallestUnique=help(transcript_id)), by="gene_id"]
##### count peptide groups
Peptides_cluster = Peptides_cluster[,.(Peptide.IDs=Peptide.IDs,transcript_id = transcript_id, smallestUnique=help(transcript_id),
                                         NuniquePeptides=sum(smallestUnique)), by= "gene_id"]

#remove mapping to more than one gene
Peptides_cluster = Peptides_cluster[!grepl(";",Peptides_cluster$gene_id),]

# decomposition of peptide list
Peptides_new = Peptides_cluster[,.(Peptide.ID=unlist(strsplit(Peptide.IDs,split=";"))), 
                                    by=c( "gene_id","Peptide.IDs","transcript_id","smallestUnique","NuniquePeptides")]

# filter intensity values from peptide ids that map to genes with at least two uniquely identifiable transcript groups
new_MS=as.data.table(merge(MS,Peptides_new,by="Peptide.ID"))
new_MS=new_MS[new_MS$NuniquePeptides>1,]
new_MS=new_MS[new_MS$smallestUnique,]

write.table(new_MS, row.names=F, col.names = T, sep = "\t", quote = F, file = out)
