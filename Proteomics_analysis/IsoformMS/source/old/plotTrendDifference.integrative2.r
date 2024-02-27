list.of.packages <- c("data.table","ggplot2","ggpubr","pheatmap","ggcorrplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("ggcorrplot"))


ratio.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/ratio.atLeastTwoProteinIsoforms.integrative2.tab"
overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/overlap.tab"
out="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/trendChanges.ONT.pdf"
nCluster=2
RNA.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/ONT-seq.txt"
args <- commandArgs(trailingOnly = TRUE) 


ratio.path = args[1]
overlap.path = args[2]
RNA.path = args[3]
nCluster=as.numeric(args[4])
method=args[5]
out =args[6]
out2=args[7]
out3=args[8]
out4=args[9]

RATIO=fread(ratio.path, header=T)
overlap=read.table(overlap.path)
RNA=fread(RNA.path, header=T)[,c("log2FoldChange","transcript_id","Comparison")]
colnames(RNA)=c("ONT_log2FC","transcript_id","Comparison")

TEMP=RATIO[,.(Protein.ID=unlist(strsplit(Protein.IDs,split=";")),
                    Gene.IDs=Gene.IDs,
              Gene.names=unique(Gene.names),
                    groupSize=length(unlist(strsplit(Protein.IDs,split=";")))),by="Protein.IDs"]

TEMP=merge(TEMP,overlap, by.x="Protein.ID", by.y="V1")
TEMP=TEMP[TEMP$V2>TEMP$groupSize,]
TEMP=TEMP[,.(Count=.N, Gene.IDs=unique(Gene.IDs)),by="Protein.IDs"]
TEMP=subset(TEMP,Count>=2)

genes_of_interest=unique(TEMP$Gene.IDs)

RATIO_FILTER=subset(RATIO, Gene.IDs %in% genes_of_interest)
############get unique gene names
RATIO_FILTER=merge(RATIO_FILTER,
                   RATIO_FILTER[,.(Gene.names.unique=paste(unique(unlist(strsplit(Gene.names,split=";"))),collapse=";")),by="Gene.IDs"],
                   by="Gene.IDs")

RATIO_FILTER=merge(RATIO_FILTER,
                        RATIO_FILTER[,.(transcript_id=unlist(strsplit(Protein.IDs,split=";"))), by="Protein.IDs"],
                        by="Protein.IDs", allow.cartesian=T)

RNA=RNA[,.(ONT_Day3_vs_Day0_ratio=ONT_log2FC[Comparison=="Day3_vs_Day0"],
       ONT_Day5_vs_Day3_ratio=ONT_log2FC[Comparison=="Day5_vs_Day3"])
,by="transcript_id"]

RATIO_FILTER=merge(RATIO_FILTER,RNA,by=c("transcript_id"))

summarize=function(X, method=method){
  if(method=="median"){
    return(median(X))
  }
  if(method=="sum"){
    return(sum(X))
  }
  if(method=="max"){
    return(X[which.max(abs(X))])
  }
}

RATIO_FILTER=RATIO_FILTER[,.(Gene.IDs=unique(Gene.IDs),
                             Gene.names.unique=unique(Gene.names.unique),
                             MS_Day3_vs_Day0=unique(Day3_vs_Day0_ratio),
                             MS_Day5_vs_Day3=unique(Day5_vs_Day3_ratio),
                             ONT_Day3_vs_Day0=summarize(ONT_Day3_vs_Day0_ratio, method=method),
                             ONT_Day5_vs_Day3=summarize(ONT_Day5_vs_Day3_ratio, method=method)
                             ),by=c("Protein.IDs")]

RATIO_FILTER_TIDY = RATIO_FILTER[,.(MS_log2FC=c(MS_Day3_vs_Day0,
                                                MS_Day5_vs_Day3
                                                ),
                                    ONT_log2FC=c(ONT_Day3_vs_Day0,
                                                      ONT_Day5_vs_Day3),
                                    Comparison=c("Day3_vs_Day0",
                                                 "Day5_vs_Day3")
                                    ),
                                 by=c("Protein.IDs","Gene.IDs","Gene.names.unique")]
#
pdf(out, width=12)
ylim=3
for(gene in genes_of_interest){
  p1=ggplot(data=subset(RATIO_FILTER_TIDY, Gene.IDs==gene), 
          aes(x=Comparison, y=MS_log2FC, fill=Protein.IDs)) +
  stat_summary(
    fun = median,
    alpha=0.6,
    size=2,
    geom = 'line',
    aes(group = Protein.IDs, colour = Protein.IDs)) +
  labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, Gene.IDs==gene)$Gene.names.unique),collapse=";")," (",gene,")",sep="")) +
  ylim(c(-1*ylim,ylim)) +
  theme(legend.title=element_blank()) +
  theme_bw()+
  theme(legend.title=element_blank(), axis.title.x=element_blank(),
        legend.position = "bottom") +
   geom_hline(yintercept=0, linetype="dashed")

  p2=ggplot(data=subset(RATIO_FILTER_TIDY, Gene.IDs==gene), 
            aes(x=Comparison, y=ONT_log2FC, fill=Protein.IDs)) +
    stat_summary(
      fun = median,
      alpha=0.6,
      size=2,
      geom = 'line',
      aes(group = Protein.IDs, colour = Protein.IDs)) +
    labs(title=paste(paste(unique(subset(RATIO_FILTER_TIDY, Gene.IDs==gene)$Gene.names.unique),collapse=";")," (",gene,")",sep="")) +
    ylim(c(-1*ylim,ylim)) +
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
}
TREND_DIFFERENCE=RATIO_FILTER[,.(MS_Day3_vs_Day0=absTrendDifference(MS_Day3_vs_Day0),
                                               MS_Day5_vs_Day3=absTrendDifference(MS_Day5_vs_Day3),
                                               ONT_Day3_vs_Day0=absTrendDifference(ONT_Day3_vs_Day0),
                                               ONT_Day5_vs_Day3=absTrendDifference(ONT_Day5_vs_Day3),
                Gene.names=paste(unique(Gene.names.unique),collapse=";")
                             ),by=c("Gene.IDs")]

TREND_DIFFERENCE$Gene.names.unique=make.unique(TREND_DIFFERENCE$Gene.names)



M=as.data.frame(TREND_DIFFERENCE)[,grepl("MS_",colnames(TREND_DIFFERENCE))]
rownames(M)=TREND_DIFFERENCE$Gene.names.unique

seed=-1
bestCluster=rep(0,times=nCluster)
while(length(unique(bestCluster))!=nCluster){
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


M_cluster=cbind(M,10*kmeans_clustering)

cluster_rows= hclust(dist(M_cluster, method="euclidean"),method="average")

annotation_row = data.frame(
  Cluster = as.factor(kmeans_clustering)
)
row.names(annotation_row)=row.names(M)

pheatmap(M,
         scale = "none",
         cluster_rows=cluster_rows,
         annotation_row = annotation_row,
         cluster_cols=F,
         show_rownames = T,
         color=rev(hcl.colors(10, "Purples 3")),
         cutree_rows = nCluster,
         treeheight_row = 0,
         filename=out2)

TREND_DIFFERENCE$cluster1=kmeans_clustering*100
TREND_DIFFERENCE=TREND_DIFFERENCE[cluster_rows$order,]
TREND_DIFFERENCE$cluster2=1:nrow(TREND_DIFFERENCE)

###############################################################################
TREND_DIFFERENCE_TIDY=TREND_DIFFERENCE[,.(MS_trend=c(MS_Day3_vs_Day0,
                                                     MS_Day5_vs_Day3),
                                          ONT_trend=c(ONT_Day3_vs_Day0,
                                                           ONT_Day5_vs_Day3),
                                          Comparison=c("Day3_vs_Day0",
                                                       "Day5_vs_Day3"),
                                          cluster1=cluster1,
                                          cluster2=cluster2),by="Gene.names.unique"]


max=round(max(abs(TREND_DIFFERENCE_TIDY$MS_trend)),0)

pdf(out3, width=12)
ggplot(data = TREND_DIFFERENCE_TIDY,aes(y=Comparison, 
                                        x = reorder(Gene.names.unique,cluster1 + cluster2),
                                        #y=reorder(name,-1*FDR), 
                                        colour=ONT_trend, 
                                        size=MS_trend)) +
  geom_point() +
  scale_fill_gradient(limits = c(0,max), breaks=seq(0,max)) +
  scale_size(range = c(2, 15)) +
  scale_colour_gradientn(colours=c("grey","red","purple")) +
  labs(x="", y="", 
       colour="Transcript Isoform",
       size="Protein Isoform") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

CORR=cor(as.data.frame(RATIO_FILTER)[,grepl("Day",colnames(RATIO_FILTER))], method = "pearson")[3:4,1:2]


pdf(out4)
ggcorrplot(CORR, method = "circle",lab = TRUE, digits=1, legend.title="Pearson\nCorrelation")
dev.off()

