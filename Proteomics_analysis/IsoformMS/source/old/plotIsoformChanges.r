list.of.packages <- c("data.table","ggplot2","ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))


MS.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3/de-novo/PeptidesMappingDistinctTranscriptGroups.txt"
gtf.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/GRCh38_p12/nanopore.SQANTI3.gtf"

overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3/de-novo/transcriptsWithUniquePeptide.overlap.tab"
minDays=2
nDataPoint=2
nTranscriptGroup=2

args <- commandArgs(trailingOnly = TRUE) 

MS.path = args[1]
gtf.path= args[2]
overlpa.path = args[3]
minDays = as.numeric(args[4])
nDataPoint =as.numeric(args[5])
nTranscriptGroup = as.numeric(args[6])
out = args[7]
out2=args[8]
#
gtf=read.table(gtf.path, sep="\t")
gtf=gtf[gtf$V3=="transcript",]
column_nine=do.call(rbind,strsplit(gtf$V9,split="; "))


index_gene_id=which(colSums(apply(column_nine,2,grepl, pattern="gene_id"))==nrow(column_nine))

index_gene_name=which(colSums(apply(column_nine,2,grepl, pattern="gene_name") | apply(column_nine,2,grepl, pattern="oId"))==nrow(column_nine))
annotation = data.table(gene_id=substring(column_nine[,index_gene_id],9),
                        gene_name=column_nine[,index_gene_name])

annotation$gene_name[grepl("gene_name",annotation$gene_name)]=substring(annotation$gene_name[grepl("gene_name",annotation$gene_name)],11)
annotation$gene_name[grepl("oId",annotation$gene_name)]=NA
annotation=unique(annotation)
annotation=annotation[,.(gene_name=paste(unique(gene_name[!is.na(gene_name)]),collapse="|")),by="gene_id"]


#prepare MS table
MS=fread(MS.path, sep="\t", header=T)

MS=merge(MS,annotation, by="gene_id")
# peptide intensities from different peptides are summarized 
# 1. sumLFQ: sum of all matching intensity values
# 2. meanLFQ: mean of all matching intensity values
new_MS=MS[,.(meanLFQ=mean(Intensity,na.rm=T), 
                medianLFQ=median(Intensity,na.rm=T), 
                sumLFQ=sum(Intensity,na.rm=T), 
                gene_name=paste(unique(gene_name),collapse="|"),
                gene_id=paste(unique(gene_id),collapse="|"),
                 Day=unique(Day)), 
              by=c("transcript_id","Experiment")]


# normalize values to 1
mean_norm_1 = new_MS[,.(Experiment=unique(Experiment), gene_id=unique(gene_id), mean_norm_1=(meanLFQ-0)/(max(meanLFQ)-0)), by="transcript_id"]
new_MS=merge(new_MS,mean_norm_1,by=c("transcript_id","gene_id","Experiment"))

sum_norm_1 = new_MS[,.(Experiment=unique(Experiment), gene_id=unique(gene_id), sum_norm_1=(sumLFQ-0)/(max(sumLFQ)-0)), by="transcript_id"]
new_MS=merge(new_MS,sum_norm_1,by=c("transcript_id","gene_id","Experiment"))

print("# isoforms")
print(length(unique(new_MS$transcript_id)))
# remove transcripts without overlap
overlap = fread(overlap.path, col.names = c("transcriptID","overlapN"))


new_T = data.table(transcriptGroup=unique(new_MS$transcript_id))
new_T=new_T[,.(temp=unlist(strsplit(transcriptGroup,split=";"))),by="transcriptGroup"]
new_T=new_T[,.(transcriptGroup=transcriptGroup,transcriptID=unlist(strsplit(temp,split="\\|"))),by="temp"]

temp=new_T[,.(groupSize=.N),by="transcriptGroup"]

select=merge(merge(new_T[,c("transcriptGroup","transcriptID")], temp, by="transcriptGroup"),overlap,by="transcriptID")

new_MS = new_MS[new_MS$transcript_id%in%unique(select$transcriptGroup[select$groupSize<select$overlapN]),]

print("# isoforms reoved non-overlapping")
print(length(unique(new_MS$transcript_id)))

# data must be available for at least minDays
select= new_MS[,.(dataP=ifelse(length(unique(Day))>=minDays,T,F)),by=c("transcript_id")]
new_MS = new_MS[new_MS$transcript_id%in%select$transcript_id[select$dataP],]

# at least nDataPoint data points must be available for each transcript, 
# at least nTranscriptGroup transcripts for each gene
select= new_MS[,.(dataP=ifelse(.N>=nDataPoint,T,F)),by=c("gene_id","transcript_id")]
select = select[,.(dataP=sum(dataP)>=nTranscriptGroup),by=c("gene_id")]
new_MS = new_MS[new_MS$gene_id%in%select$gene_id[select$dataP],]

print("# isoforms after data point anf transcript number filter")
print(length(unique(new_MS$transcript_id)))

write.table(new_MS, row.names=F, col.names=T, quote = F, sep = "\t", file=out)

pdf(out2,width=12)
for (gene in unique(new_MS$gene_id)){
  selection=new_MS[new_MS$gene_id==gene,]
  gene_name=unique(selection$gene_name)
  
  p1=ggplot(data=selection, aes(x=Day, y=sum_norm_1, fill =transcript_id)) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id),
      position = position_dodge(width = 0.9) #this has to be added
    ) +
    geom_boxplot(position=position_dodge(1), width = 0.5) +
    ylim(0, 1) +
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize=0.5) +
    labs(title=paste(gene_name," (",gene,")",sep="")) +
    theme(legend.title=element_blank()) +
    theme_bw()
  
  p2=ggplot(data=selection, aes(x=Day, y= mean_norm_1, fill =transcript_id)) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id),
      position = position_dodge(width = 0.9) #this has to be added
    ) +
    geom_boxplot(position=position_dodge(1), width = 0.5) +
    ylim(0, 1) +
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize=0.5) +
    labs(title=paste(gene_name," (",gene,")",sep="")) +
    theme(legend.title=element_blank()) +
    theme_bw()
  
  
  p=ggarrange(p1,p2, 
              nrow = 1, ncol=2,
              align="hv",
              common.legend = TRUE, legend="top")
  
  print(p)
}
dev.off()

