list.of.packages <- c("data.table","ggplot2","ggpubr","pheatmap","ggcorrplot","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("ggcorrplot"))
suppressPackageStartupMessages(library("tidyr"))


ratio.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/MSratio.atLeastTwoProteinIsoforms.Illumina.tab"
overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/overlap.tab"
lfq.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/FlashLFQ/FashLFQ.txt"
experimentalDesign.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/ExperimentalDesign.continuous.txt"

nCluster=2
RNA.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/RNA-seq.txt"
out.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/trendChanges.Illumina."
log2FC_threshold=1
args <- commandArgs(trailingOnly = TRUE) 


ratio.path = args[1]
overlap.path = args[2]
RNA.path = args[3]
nCluster=as.numeric(args[4])
out.path =args[5]

out=paste(out.path,"png",sep="")
out2=paste(out.path,"heatmap.pdf",sep="")
out3=paste(out.path,"bubblePlot.pdf",sep="")
out4=paste(out.path,"corr.all.pdf",sep="")
out5=paste(out.path,"corr.trendDifference.pdf",sep="")

RATIO=fread(ratio.path, header=T)
overlap=read.table(overlap.path)
RNA=fread(RNA.path, header=T)[,c("log2FoldChange","transcript_id","Comparison")]
colnames(RNA)=c("RNA_log2FC","transcript_id","Comparison")

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
RATIO_FILTER_TIDY=merge(RATIO_FILTER_TIDY,RNA,by=c("transcript_id","Comparison"))

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

#########



png(out) #, width=12, units="in")

for(gene in genes_of_interest){
  ylim_MS=ceiling(max(abs(subset(RATIO_FILTER_TIDY, gene_id==gene)$MS_log2FC)))
  ylim_RNA=ceiling(max(abs(subset(RATIO_FILTER_TIDY, gene_id==gene)$RNA_log2FC)))
  
  name=unique(subset(RATIO_FILTER_TIDY, gene_id==gene)$gene_name)
  p1=ggplot(data=subset(RATIO_FILTER_TIDY, gene_id==gene), 
          aes(x=Comparison, y=MS_log2FC, fill=transcript_id)) +
  stat_summary(
    fun = median,
    alpha=0.6,
    size=2,
    geom = 'line',
    aes(group = transcript_id, colour = transcript_id)) +
  labs(title=paste(paste(name,collapse=";")," (",gene,")",sep="")) +
  ylim(c(-1*ylim_MS,ylim_MS)) +
  theme(legend.title=element_blank()) +
  theme_bw()+
  theme(legend.title=element_blank(), axis.title.x=element_blank(),
        legend.position = "bottom") +
   geom_hline(yintercept=0, linetype="dashed")

  p2=ggplot(data=subset(RATIO_FILTER_TIDY, gene_id==gene), 
            aes(x=Comparison, y=RNA_log2FC, fill=transcript_id)) +
    stat_summary(
      fun = median,
      alpha=0.6,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)) +
    labs(title=paste(paste(name,collapse=";")," (",gene,")",sep="")) +
    ylim(c(-1*ylim_RNA,ylim_RNA)) +
    theme(legend.title=element_blank()) +
    theme_bw()+
    theme(legend.title=element_blank(), axis.title.x=element_blank(),
          legend.position = "bottom") +
    geom_hline(yintercept=0, linetype="dashed")
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

TREND_DIFFERENCE_TIDY=as.data.table(RATIO_FILTER_TIDY)[,.(MS_trend=absTrendDifference(MS_log2FC),
                                                          RNA_trend=absTrendDifference(RNA_log2FC),
                                                          Gene.names=make.unique(rep(paste(unique(gene_name),collapse=";"),
                                                                                     times=length(absTrendDifference(MS_log2FC))))
),by=c("gene_id","Comparison")]

myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

TREND_DIFFERENCE=myspread(TREND_DIFFERENCE_TIDY,Comparison, c(MS_trend,RNA_trend))

row.names(TREND_DIFFERENCE)=TREND_DIFFERENCE$Gene.names
genes_highTrendChange=unique(TREND_DIFFERENCE_TIDY$gene_id[TREND_DIFFERENCE_TIDY$MS_trend>log2FC_threshold])

M=as.data.frame(TREND_DIFFERENCE)[TREND_DIFFERENCE$gene_id %in% genes_highTrendChange,grepl("MS_",colnames(TREND_DIFFERENCE))]

seed=-1
bestCluster=rep(0,times=ncol(M))
while(length(unique(bestCluster))<min(ncol(M),nCluster)){
  seed=seed+1
  print(seed)
  set.seed(seed)
  kmeans = pheatmap(M,
                    scale = "row",
                    kmeans_k=nCluster,
                    cluster_cols=F,
                    cluster_rows=F
  )
  bestCluster=apply(kmeans$kmeans$centers,2,which.max)
}

help=function(X,bestCluster=bestCluster){
  return(which(X==bestCluster))
}
if(length(unique(bestCluster))==length(bestCluster)){
  kmeans_clustering=unlist(lapply(kmeans$kmeans$cluster, help, bestCluster=unname(bestCluster)))
}else{
  kmeans_clustering=kmeans$kmeans$cluster
  
}


M_cluster=cbind(M,cluster_1=100*kmeans_clustering)

cluster_rows= hclust(dist(M_cluster, method="euclidean"),method="average")

M$cluster1=kmeans_clustering*100
M=M[cluster_rows$order,]
M$cluster2=1:nrow(M)
M=M[order(M$cluster1+M$cluster2),]
M$order=1:nrow(M)

annotation_row = data.frame(
  Cluster = as.factor(sort(kmeans_clustering))
)
row.names(annotation_row)=row.names(M)

pheatmap(M[,grepl("MS_trend",colnames(M))],
         scale = "none",
         cluster_rows=F,
         annotation_row = annotation_row,
         cluster_cols=F,
         show_rownames = T,
         color=rev(hcl.colors(10, "Purples 3")),
         cutree_rows = nCluster,
         treeheight_row = 0,
         filename = out2)

M$Gene.names=row.names(M)

###############################################################################
TREND_DIFFERENCE_TIDY=merge(TREND_DIFFERENCE_TIDY,M[,c("Gene.names","order")], by="Gene.names")


max=round(max(abs(TREND_DIFFERENCE_TIDY$MS_trend)),0)

pdf(out3, width=12)
ggplot(data = TREND_DIFFERENCE_TIDY,aes(y=Comparison, 
                                        x = reorder(Gene.names,order),
                                        #y=reorder(name,-1*FDR), 
                                        colour=RNA_trend, 
                                        size=MS_trend)) +
  geom_point() +
  scale_fill_gradient(limits = c(0,max), breaks=seq(0,max)) +
  scale_size(range = c(2, 10)) +
  scale_colour_gradient(low="grey",high="red") +
  labs(x="", y="", 
       colour="Transcript Isoform",
       size="Protein Isoform") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

RATIO_FILTER=myspread(RATIO_FILTER_TIDY,Comparison, c(MS_log2FC,RNA_log2FC))
colnames(RATIO_FILTER)=gsub("_log2FC","",colnames(RATIO_FILTER))
RATIO_FILTER=RATIO_FILTER[ , c(grep("MS",colnames(RATIO_FILTER)),grep("RNA",colnames(RATIO_FILTER)))]
CORR=cor(as.data.frame(RATIO_FILTER)[,grepl("Day",colnames(RATIO_FILTER))], method = "spearman")

CORR=CORR[grep("RNA",row.names(CORR)),grep("MS",row.names(CORR))]
pdf(out4)
ggcorrplot(CORR, method = "circle",lab = TRUE, digits=1, legend.title="Spearman\nCorrelation")
dev.off()

RATIO_FILTER=myspread(RATIO_FILTER_TIDY,Comparison, c(MS_log2FC,RNA_log2FC))
RATIO_FILTER=subset(RATIO_FILTER, gene_id %in% genes_highTrendChange)
colnames(RATIO_FILTER)=gsub("_log2FC","",colnames(RATIO_FILTER))
RATIO_FILTER=RATIO_FILTER[ , c(grep("MS",colnames(RATIO_FILTER)),grep("RNA",colnames(RATIO_FILTER)))]
CORR=cor(as.data.frame(RATIO_FILTER)[,grepl("Day",colnames(RATIO_FILTER))], method = "spearman")

CORR=CORR[grep("RNA",row.names(CORR)),grep("MS",row.names(CORR))]
pdf(out5)
ggcorrplot(CORR, method = "circle",lab = TRUE, digits=1, legend.title="Spearman\nCorrelation")
dev.off()


