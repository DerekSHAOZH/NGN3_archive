list.of.packages <- c("data.table","ggplot2","ggpubr","pheatmap","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("tidyr"))


ratio.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/MSratio.atLeastTwoProteinIsoforms.tab"

lfq.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/FlashLFQ/FashLFQ.txt"
experimentalDesign.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/ExperimentalDesign.continuous.txt"

overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/overlap.tab"
nCluster=5

log2FC_threshold=1
args <- commandArgs(trailingOnly = TRUE) 


ratio.path = args[1]
lfq.path=args[2]
overlap.path = args[3]
nCluster=as.numeric(args[4])
log2FC_threshold=as.numeric(args[5])
experimentalDesign.path=args[6]
out = args[7]
out2 =args[8]

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
RATIO_FILTER_TIDY=gather(RATIO_FILTER, key=Comparison, value=log2FC, colnames(RATIO_FILTER)[grepl("Day",colnames(RATIO_FILTER))])

#write.table(unique(RATIO_FILTER$gene_name), col.names = F, row.names = F,quote=F,file="/home/bressin/Downloads/list.txt")
####################
#SD=read.table(sd.path, header=T)
#SD=subset(SD, transcript_id %in% transcripts_of_interest)
#SD_TIDY=gather(SD, key=Comparison, value=sd, colnames(SD)[grepl("Day",colnames(SD))])

########################
#RATIO_FILTER_TIDY=merge(RATIO_FILTER_TIDY,SD_TIDY, bx=c("transcript_id","gene_id","gene_name","Comparison"))
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
               unique(RATIO_FILTER_TIDY[,c("transcript_id","gene_id","gene_name")]),
               by.x="Protein.Groups", by.y="transcript_id")

min_max=function(X,min=0){return((X-min)/(max(X)-min))}

LFQ_TIDY=LFQ_TIDY[,.(Condition=Condition,
            FlashLFQ=min_max(FlashLFQ),
            Biorep=Biorep,
            gene_id=gene_id,
            gene_name=gene_name),by="Protein.Groups"]
####




pdf(out,width=12)
for(gene in genes_of_interest){
  
  ylim=ceiling(max(abs(subset(RATIO_FILTER_TIDY, gene_id==gene)$log2FC)))
  
  
  
  p1=ggplot(data=subset(LFQ_TIDY, gene_id==gene), 
            aes(x=as.factor(Condition), y=FlashLFQ, fill=Protein.Groups)) +
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize=0.5) +
    labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name),collapse=";")," (",gene,")",sep="")) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = Protein.Groups, colour = Protein.Groups),
      position = position_dodge(width = 0.9) #this has to be added
    ) +
    geom_boxplot(position=position_dodge(1), width = 0.5, alpha=0.5) +
    theme_bw()+
    ylab("LFQ (a.u.)") +
    theme(legend.title=element_blank(), axis.title.x=element_blank(),
          legend.position = "bottom") 
  
  
  p2=ggplot(data=subset(RATIO_FILTER_TIDY, gene_id==gene), 
            aes(x=Comparison, y=log2FC, fill=transcript_id)) +
    stat_summary(
      fun = median,
      alpha=0.6,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)) +
    labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name),collapse=";")," (",gene,")",sep="")) +
    ylim(c(-1*ylim,ylim)) +
    theme(legend.title=element_blank()) +
    theme_bw()+
    theme(legend.title=element_blank(), axis.title.x=element_blank(),
          legend.position = "bottom") +
    geom_hline(yintercept=0, linetype="dashed") #+
#    geom_errorbar(aes(ymin=log2FC-sd, ymax=log2FC+sd), width=.2)

  
  p=ggarrange(p1,p2, 
              nrow = 1, ncol=2,
              align="hv",
              common.legend = TRUE, legend="top")
  
  print(p)
}
dev.off()
  
absTrendDifference=function(X){
  if(length(X)==2){
    return(abs(X[1]-X[2]))
  }
  if(length(X)==3){
    return(c(abs(X[1]-X[2]),abs(X[1]-X[3]),abs(X[2]-X[3])))
  }
  if(length(X)==4){
    return(c(abs(X[1]-X[2]),abs(X[1]-X[3]),abs(X[1]-X[4]),
             abs(X[2]-X[3]),abs(X[2]-X[4]),abs(X[3]-X[4])))
  }
}

getTranscriptNames=function(X){
  if(length(X)==2){
    return(paste(X,collapse=";"))
  }
  if(length(X)==3){
    return(c(paste(X[1],X[2],sep=";"),paste(X[1],X[3],sep=";"),paste(X[2],X[3],sep=";")))
  }
  if(length(X)==4){
    return(c(paste(X[1],X[2],sep=";"),paste(X[1],X[3],sep=";"),paste(X[1],X[4],sep=";"),
             paste(X[2],X[3],sep=";"),paste(X[2],X[4],sep=";"),paste(X[3],X[4],sep=";")))
  }
}
TREND_DIFFERENCE_TIDY=as.data.table(RATIO_FILTER_TIDY)[,.(
                                          transcript_ids=getTranscriptNames(transcript_id),
                                          MS_trend=absTrendDifference(log2FC),
                                          Gene.names=make.unique(rep(paste(unique(gene_name),collapse=";"),
                                                                      times=length(absTrendDifference(log2FC))))
),by=c("gene_id","Comparison")]

TREND_DIFFERENCE=spread(TREND_DIFFERENCE_TIDY, key=Comparison,value=MS_trend)

genes_highTrendChange=unique(TREND_DIFFERENCE_TIDY$gene_id[round(TREND_DIFFERENCE_TIDY$MS_trend,1)>=log2FC_threshold])

transcripts_highTrendChange=unique(unlist(strsplit(TREND_DIFFERENCE_TIDY$transcript_ids[round(TREND_DIFFERENCE_TIDY$MS_trend,1)>=log2FC_threshold],split=";")))

######################### stat
print("number of genes with at least two detectable protein isoforms")
length(genes_of_interest)
print("number of genes with at high trend difference")
length(genes_highTrendChange)

################
row.names(TREND_DIFFERENCE)=TREND_DIFFERENCE$Gene.names

M=TREND_DIFFERENCE[TREND_DIFFERENCE$gene_id %in% genes_highTrendChange,
                   grepl("Day",colnames(TREND_DIFFERENCE))]
print("M rows")
nrow(M)
z_score=function(X){  return((X-mean(X))/sd(X))}
#MATRIX=t(apply(M,1,z_score))
MATRIX=M

seed=0
bestCluster=rep(0,times=nCluster)
while(length(unique(bestCluster))<min(ncol(M),nCluster)){
  seed=seed+1
  set.seed(seed)
  
  kmeans = pheatmap(MATRIX,
                    scale="row",
                    kmeans_k=nCluster,
                    cluster_cols=F,
                    cluster_rows=F
  )
  
  bestCluster=apply(kmeans$kmeans$centers,2,which.max)
  
}

print(seed)

help=function(X,bestCluster=bestCluster){
  return(which(X==bestCluster))
}
if(length(unique(bestCluster))==length(bestCluster)){
  kmeans_clustering=unlist(lapply(kmeans$kmeans$cluster, help, bestCluster=unname(bestCluster)))
}else{
  kmeans_clustering=kmeans$kmeans$cluster
  
}

MATRIX_cluster=cbind(M,cluster_1=100*kmeans_clustering)

cluster_rows= hclust(dist(MATRIX_cluster, method="euclidean"), method="average")

MATRIX_cluster=cbind(MATRIX_cluster,cluster_2=cluster_rows$order)

M_order=M[order(rowSums(MATRIX_cluster[,c("cluster_1","cluster_2")]),decreasing=F),]
annotation_row = data.frame(
  Cluster = as.factor(sort(kmeans_clustering,decreasing = F))
)
row.names(annotation_row)=row.names(M_order)

breaks=seq(0,3,0.1)

pheatmap(M_order,
         breaks=breaks,
         scale = "none",
         annotation_row = annotation_row,
         cluster_rows=F,
         cluster_cols=F,
         show_rownames = T,
         color=colorRampPalette(rev(hcl.colors(10, "Purples 3")))(length(breaks)+1),
         gaps_row = cumsum(table(kmeans_clustering)),
         treeheight_row = 0,
         widthInch=12,
         fontsize_row=12,
         fontsize_col=12,
         cellwidth=10,
        filename = out2)

###############################################################################


