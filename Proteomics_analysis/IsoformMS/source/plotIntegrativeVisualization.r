list.of.packages <- c("data.table","ggplot2","ggpubr","pheatmap","ggcorrplot","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("ggcorrplot"))
suppressPackageStartupMessages(library("tidyr"))


ratio.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/MSratio.atLeastTwoProteinIsoforms.tab"

lfq.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/FlashLFQ/FashLFQ.txt"
overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/overlap.tab"
experimentalDesign.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/ExperimentalDesign.continuous.txt"

rna.identity.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/RNA-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt"
ont.identity.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/ONT-seq_NeuroDifferentation_NGN3/isoform/countTable.RLE.txt"

out="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/integrative.pdf"

args <- commandArgs(trailingOnly = TRUE) 


ratio.path = args[1]
lfq.path=args[2]
overlap.path = args[3]
experimentalDesign.path = args[4]
rna.identity.path=args[5]
ont.identity.path=args[6]
out=args[7]

RATIO=fread(ratio.path, header=T)
overlap=read.table(overlap.path)


TEMP=RATIO[,.(transcript_id_split=unlist(strsplit(transcript_id,split=";")),
              gene_id=gene_id,
              gene_name=unique(gene_name),
              groupSize=length(unlist(strsplit(transcript_id,split=";")))),by="transcript_id"]

TEMP=merge(TEMP,overlap, by.x="transcript_id_split", by.y="V1")
TEMP=TEMP[TEMP$V2>TEMP$groupSize,]
TEMP=merge(TEMP,TEMP[,.(Count=.N),by="gene_id"], by="gene_id")
TEMP=subset(TEMP,Count>=2)

genes_of_interest=unique(TEMP$gene_id)
transcripts_of_interest=unique(TEMP$transcript_id)

RATIO_FILTER=subset(RATIO, transcript_id %in% transcripts_of_interest)

####
RATIO_FILTER_TIDY=gather(RATIO_FILTER, key=Comparison, value=MS_log2FC, colnames(RATIO_FILTER)[grepl("Day",colnames(RATIO_FILTER))])
annotation=unique(RATIO_FILTER_TIDY[,c("transcript_id","gene_id","gene_name")])
####
LFQ=read.table(lfq.path,sep="\t",header=T)
experimentalDesign=read.table(experimentalDesign.path,sep="\t", header=T)

LFQ=subset(LFQ, Protein.Groups %in% transcripts_of_interest)
LFQ=LFQ[,grepl("Intensity_",colnames(LFQ)) | grepl("Protein.Groups",colnames(LFQ))]
LFQ_TIDY=gather(LFQ, key=FileName, value=FlashLFQ, colnames(LFQ)[grepl("Intensity_",colnames(LFQ))])

LFQ_TIDY$FileName=gsub("\\.","-",substring(LFQ_TIDY$FileName,11))
LFQ_TIDY=merge(LFQ_TIDY, experimentalDesign[,c("FileName","Condition","Biorep","Fraction")], by="FileName")

LFQ_TIDY=as.data.table(LFQ_TIDY)[,.(FlashLFQ=sum(FlashLFQ)),by=c("Protein.Groups","Condition","Biorep")]
LFQ_TIDY=merge(LFQ_TIDY, 
               annotation,
               by.x="Protein.Groups", by.y="transcript_id")

min_max=function(X,min=0){return((X-min)/(max(X)-min))}

LFQ_TIDY=LFQ_TIDY[,.(Condition=Condition,
                     FlashLFQ=min_max(FlashLFQ),
                     Biorep=Biorep,
                     gene_id=gene_id,
                     gene_name=gene_name),by="Protein.Groups"]

#########

RNA_identity=read.table(rna.identity.path, header=T,sep=" ")
RNA_identity_TIDY=gather(RNA_identity, key=Sample, value=RLE, colnames(RNA_identity)[grepl("Day",colnames(RNA_identity))])

RNA_identity_TIDY$Day=do.call(rbind,strsplit(RNA_identity_TIDY$Sample,split="_"))[,1]
RNA_identity_TIDY$Replicate=do.call(rbind,strsplit(RNA_identity_TIDY$Sample,split="_"))[,2]
RNA_identity_TIDY=merge(annotation,
                        RNA_identity_TIDY,  by.y="X", by.x="transcript_id")

RNA_identity_TIDY=as.data.table(RNA_identity_TIDY)[,.(Day=Day,
                     RLE=min_max(RLE),
                     Replicate=Replicate,
                     gene_id=gene_id,
                     gene_name=gene_name),by="transcript_id"]

#########

ONT_identity=read.table(ont.identity.path, header=T,sep=" ")
ONT_identity_TIDY=gather(ONT_identity, key=Sample, value=RLE, colnames(ONT_identity)[grepl("Day",colnames(ONT_identity))])

ONT_identity_TIDY$Day=do.call(rbind,strsplit(ONT_identity_TIDY$Sample,split="_"))[,1]
ONT_identity_TIDY$Replicate=do.call(rbind,strsplit(ONT_identity_TIDY$Sample,split="_"))[,2]
ONT_identity_TIDY=merge(annotation,
                        ONT_identity_TIDY,  by.y="X", by.x="transcript_id")

ONT_identity_TIDY=as.data.table(ONT_identity_TIDY)[,.(Day=Day,
                                                      RLE=min_max(RLE),
                                                      Replicate=Replicate,
                                                      gene_id=gene_id,
                                                      gene_name=gene_name),by="transcript_id"]
########################

pdf(out, width=10, height = 4)

for(gene in genes_of_interest){
  col=c("#F8766D","#00BFC4","#619CFF")
  name=unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name)
  p1=ggplot(data=subset(LFQ_TIDY, gene_id==gene & Protein.Groups!="TCONS_00001947"), 
            aes(x=as.factor(Condition), y=FlashLFQ, fill=Protein.Groups)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + #position=position_dodge(1),
    labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name),collapse=";")," (",gene,")",sep="")) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = Protein.Groups, colour = Protein.Groups)) + 
    #  position = position_dodge(width = 0.9) #this has to be added
    #) +
    geom_boxplot(width = 0.5, alpha=0.5) + #position=position_dodge(1), 
    ylim(0,1) +
    theme_bw()+
    ylab("LFQ (a.u.)") +
    xlab("Protein Isoform (MS)") +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    theme(legend.title=element_blank(),
          legend.position = "bottom") 
  
  p2=ggplot(data=subset(RNA_identity_TIDY, gene_id==gene & transcript_id!="TCONS_00001947"), 
            aes(x=as.factor(Day), y=RLE, fill=transcript_id)) +
    geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5) + #position=position_dodge(1), 
    labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name),collapse=";")," (",gene,")",sep="")) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)) +
     # position = position_dodge(width = 0.9) #this has to be added
    #) +
    geom_boxplot(width = 0.5, alpha=0.5) + #
    ylim(0,1) +
    theme_bw()+
    ylab("RLE (a.u.)") +
    xlab("Gene Transcript Isoform (Illumina)") +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    theme(legend.title=element_blank(), 
          legend.position = "bottom") 

  p3=ggplot(data=subset(ONT_identity_TIDY, gene_id==gene & transcript_id!="TCONS_00001947"), 
            aes(x=as.factor(Day), y=RLE, fill=transcript_id)) +
    geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5) + #position=position_dodge(1), 
    labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name),collapse=";")," (",gene,")",sep="")) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)) +
     # position = position_dodge(width = 0.9) #this has to be added
    #) +
    geom_boxplot(width = 0.5, alpha=0.5) + #position=position_dodge(1), 
    ylim(0,1) +
    theme_bw()+
    ylab("RLE (a.u.)") +
    xlab("Gene Transcript Isoform (ONT)") +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    theme(legend.title=element_blank(),
          legend.position = "bottom") 
  
  p=ggarrange(p1,p2,p3, 
              nrow = 1, ncol=3,
              align="hv",
              common.legend = TRUE, legend="top")
  
  print(p)
}
dev.off()
