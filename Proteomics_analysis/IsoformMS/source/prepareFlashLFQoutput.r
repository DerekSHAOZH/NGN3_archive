list.of.packages <- c("data.table","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("tidyr"))

data.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/FlashLFQ/BayesianFoldChangeAnalysis.all.tsv"
annotation.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/reference_day0/annotation.info.txt"
mode="atLeastTwoProteinIsoforms"
nCol=14
compareL="Day1_vs_Day0;Day2_vs_Day1;Day3_vs_Day2;Day4_vs_Day3;Day5_vs_Day4"
args <- commandArgs(trailingOnly = TRUE) 

data.path = args[1]
annotation.path = args[2]
nCol = as.numeric(args[3])
mode =args[4]
out = args[5]
out2 = args[6]
compareL=args[7]

compare=unlist(strsplit(compareL,split=";"))

data=read.table(data.path,sep="\t", header = T)
#remove seperation column X
data=data[,!grepl("X",colnames(data))]

data=split.default(data,rep(1:(ncol(data)/nCol),each=nCol))
data=rbindlist(data, use.names=F)

data=subset(data, !grepl("CON__",Protein.Group))
data=subset(data, Number.of.Peptides>0)

data$Condition=paste(data$Treatment.Condition,data$Control.Condition,sep="_vs_")

# select only comparisons of interest
data=subset(data,Condition %in% compare)


FDR=spread(data[,c("Protein.Group","False.Discovery.Rate","Condition")], key=Condition,value=False.Discovery.Rate)
log2FC=spread(data[,c("Protein.Group","Protein.Log2.Fold.Change","Condition")], key=Condition,value=Protein.Log2.Fold.Change)

#add gene id

annotation=read.table(annotation.path, header =T)
log2FC=merge(annotation,log2FC, by.y="Protein.Group",by.x="transcript_id")
FDR=merge(annotation,FDR, by.y="Protein.Group",by.x="transcript_id")

genes_of_interest=names(table(log2FC$gene_id)[table(log2FC$gene_id)>=2])

if(mode=="atLeastTwoProteinIsoforms"){
  log2FC=subset(log2FC,gene_id %in% genes_of_interest)
  FDR=subset(FDR,gene_id %in% genes_of_interest)
}

write.table(log2FC, row.names = F, col.names = T, sep = "\t", quote = F, file=out)
write.table(FDR, row.names = F, col.names = T, sep = "\t", quote = F, file=out2)

