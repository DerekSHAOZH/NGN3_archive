list.of.packages <- c("data.table","ggplot2","ggpubr","pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("pheatmap"))


prefix="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/de-novo/LFQ.all."
MS.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/de-novo/LFQ.all.tab"
overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/de-novo/all.overlap.tab"


summarizeLFQMethod="mean"
args <- commandArgs(trailingOnly = TRUE) 


MS.path = args[1]
overlap.path = args[2]
prefix = args[3]
summarizeLFQMethod = args[4]

TABLE=as.data.table(read.table(MS.path, header=1))
overlap=read.table(overlap.path)

TABLE_TEMP=TABLE[,.(Protein.ID=unlist(strsplit(Protein.IDs,split=";")),
                    Gene.IDs=Gene.IDs,
                    groupSize=length(unlist(strsplit(Protein.IDs,split=";")))),by="Protein.IDs"]

TABLE_TEMP=merge(TABLE_TEMP,overlap, by.x="Protein.ID", by.y="V1")
TABLE_TEMP=TABLE_TEMP[TABLE_TEMP$V2>TABLE_TEMP$groupSize,]
TABLE_TEMP=TABLE_TEMP[,.(Count=.N, Gene.IDs=unique(Gene.IDs)),by="Protein.IDs"]
TABLE_TEMP=subset(TABLE_TEMP,Count>=2)

gene_of_interest=unique(TABLE_TEMP$Gene.IDs)

TABLE_TIDY=TABLE[,.(LFQ=c(LFQ.intensity.Day0.B,LFQ.intensity.Day0.C,LFQ.intensity.Day0.D,
                         LFQ.intensity.Day1.B,LFQ.intensity.Day1.C,LFQ.intensity.Day1.D,
                         LFQ.intensity.Day2.B,LFQ.intensity.Day2.C,LFQ.intensity.Day2.D,
                         LFQ.intensity.Day3.B,LFQ.intensity.Day3.C,LFQ.intensity.Day3.D,              
                         LFQ.intensity.Day4.B,LFQ.intensity.Day4.C,LFQ.intensity.Day4.D,              
                         LFQ.intensity.Day5.B.C,LFQ.intensity.Day5.D.E,LFQ.intensity.Day5.F.G),
                   Replicate=rep(1:3,times=6),
                   Day=rep(paste("Day",0:5),each=3),
                   Gene.IDs=Gene.IDs,
                   Gene.NAMEs=Gene.NAMEs),by="Protein.IDs"]

# remove non-overlapping transcripts
TABLE_TIDY=TABLE_TIDY[TABLE_TIDY$Gene.IDs %in% gene_of_interest,]

#remove lowly expressed thrid, fourth transcripts 
TEMP=TABLE_TIDY[,.(meanLFQ=mean(LFQ)),by=c("Gene.IDs","Protein.IDs")]
TEMP=merge(TEMP,TEMP[,.(count=.N),by="Gene.IDs"],by="Gene.IDs")

TEMP=merge(TEMP,TEMP[,.(Protein.IDs=Protein.IDs[order(meanLFQ, decreasing = T)],
                        ExpressinOrder=1:.N),by="Gene.IDs"],by="Protein.IDs")

removeTranscript=TEMP$Protein.IDs[TEMP$count>2 & TEMP$ExpressinOrder>2]

###########################
TABLE_TIDY=TABLE_TIDY[!(TABLE_TIDY$Protein.IDs %in% removeTranscript),]

computeLog2FC = function(x,y){
  pc=1
  return(log2((x+pc)/(y+pc)))
}

summarizeLFQ=function(X, summarizeLFQMethod=summarizeLFQMethod){
  if(summarizeLFQMethod=="mean"){
    return(mean(X))
  }
  
  if(summarizeLFQMethod=="median"){
    return(median(X))
  }
  
}
TABLE_TREND_TIDY=TABLE_TIDY[,.(LFQ=summarizeLFQ(LFQ, summarizeLFQMethod=summarizeLFQMethod),
                               Gene.NAMEs=unique(Gene.NAMEs),
                               Gene.IDs=unique(Gene.IDs)),by=c("Protein.IDs","Day")]

TABLE_TREND_TIDY=TABLE_TIDY[,.(t_test=c(t.test(LFQ[Day=="Day 1"], LFQ[Day=="Day 0"], alternative = "two.sided")$p.value,
                                      t.test(LFQ[Day=="Day 2"], LFQ[Day=="Day 1"], alternative = "two.sided")$p.value,
                                      t.test(LFQ[Day=="Day 3"], LFQ[Day=="Day 2"], alternative = "two.sided")$p.value,
                                      t.test(LFQ[Day=="Day 4"], LFQ[Day=="Day 3"], alternative = "two.sided")$p.value,
                                      t.test(LFQ[Day=="Day 5"], LFQ[Day=="Day 4"], alternative = "two.sided")$p.value),
                             log2FC=c(computeLog2FC(summarizeLFQ(LFQ[Day=="Day 1"], summarizeLFQMethod=summarizeLFQMethod),summarizeLFQ(LFQ[Day=="Day 0"], summarizeLFQMethod=summarizeLFQMethod)),
                                      computeLog2FC(summarizeLFQ(LFQ[Day=="Day 2"], summarizeLFQMethod=summarizeLFQMethod),summarizeLFQ(LFQ[Day=="Day 1"], summarizeLFQMethod=summarizeLFQMethod)),
                                      computeLog2FC(summarizeLFQ(LFQ[Day=="Day 3"], summarizeLFQMethod=summarizeLFQMethod),summarizeLFQ(LFQ[Day=="Day 2"], summarizeLFQMethod=summarizeLFQMethod)),
                                      computeLog2FC(summarizeLFQ(LFQ[Day=="Day 4"], summarizeLFQMethod=summarizeLFQMethod),summarizeLFQ(LFQ[Day=="Day 3"], summarizeLFQMethod=summarizeLFQMethod)),
                                      computeLog2FC(summarizeLFQ(LFQ[Day=="Day 5"], summarizeLFQMethod=summarizeLFQMethod),summarizeLFQ(LFQ[Day=="Day 4"], summarizeLFQMethod=summarizeLFQMethod))),
                             Comparison=c("Day 1 vs 0","Day 2 vs 1","Day 3 vs 2","Day 4 vs 3","Day 5 vs 4"),
                             Gene.NAMEs=unique(Gene.NAMEs),
                             Gene.IDs=unique(Gene.IDs)),
                          by=c("Protein.IDs")]

TABLE_TREND_TIDY=merge(TABLE_TREND_TIDY,
                       TABLE_TREND_TIDY[,.(absTrendDifference=abs(log2FC[1]-log2FC[2])),by=c("Gene.IDs","Comparison")],
                       by=c("Gene.IDs","Comparison"))

# select only examples where:
# protein expression changed significantly for at least one protein isoform
# log2 trend difference of 2 or more
pValue_thr=0.05
log2FC_thr=1

atLeastOneSigDifference=unique(TABLE_TREND_TIDY$Gene.NAMEs[TABLE_TREND_TIDY$t_test<pValue_thr &
                                                           !is.na(TABLE_TREND_TIDY$t_test)])
trendDifference=unique(TABLE_TREND_TIDY$Gene.NAMEs[TABLE_TREND_TIDY$absTrendDifference>=log2FC_thr])
atLeastOneSigDifference_trendDifference=unique(TABLE_TREND_TIDY$Gene.NAMEs[TABLE_TREND_TIDY$t_test<pValue_thr & 
                                                                           !is.na(TABLE_TREND_TIDY$t_test) & 
                                                                           TABLE_TREND_TIDY$absTrendDifference>=log2FC_thr])

length(atLeastOneSigDifference)
length(trendDifference)
length(atLeastOneSigDifference_trendDifference)

pdf(paste(prefix,"pdf",sep=""), width=12)
for(gene in gene_of_interest){

  p1=ggplot(data=subset(TABLE_TIDY, Gene.IDs==gene), 
           aes(x=Day, y=LFQ, fill=Protein.IDs)) +
    stat_summary(
      fun = median,
      size=2,
      geom = 'line',
      alpha=0.6,
      aes(group = Protein.IDs, colour = Protein.IDs),
      position = position_dodge(width = 0.9) #this has to be added
    ) +
    geom_boxplot(position=position_dodge(1), width = 0.5) +
    geom_dotplot(binaxis='y', stackdir='center',position=position_dodge(1), dotsize=0.5) +
    labs(title=paste(paste(unique(subset(TABLE_TIDY, Gene.IDs==gene)$Gene.NAMEs),collapse=";")," (",gene,")",sep="")) +
    theme(legend.title=element_blank()) +
    theme_bw()+
    theme(legend.title=element_blank(), axis.title.x=element_blank(),
          legend.position = "bottom")
  
  p2=ggplot(data=subset(TABLE_TREND_TIDY, Gene.IDs==gene), 
            aes(x=Comparison, y=log2FC, fill=Protein.IDs)) +
    stat_summary(
      fun = median,
      alpha=0.6,
      size=2,
      geom = 'line',
      aes(group = Protein.IDs, colour = Protein.IDs)) +
    labs(title=paste(paste(unique(subset(TABLE_TREND_TIDY, Gene.IDs==gene)$Gene.NAMEs),collapse=";")," (",gene,")",sep="")) +
    ylim(ifelse(max(abs(subset(TABLE_TREND_TIDY, Gene.IDs==gene)$log2FC))>2.5,25,2.5)*c(-1,1)) +
    theme(legend.title=element_blank()) +
    theme_bw()+
    theme(legend.title=element_blank(), axis.title.x=element_blank(),
          legend.position = "bottom")
  
  p=ggarrange(p1,p2, 
              nrow = 1, ncol=2,
              align="hv",
              common.legend = TRUE, legend="top")
  
  print(p)
}
dev.off()

TABLE_MATRIX=as.data.frame(TABLE_TREND_TIDY[,.(Gene.NAMEs=unique(Gene.NAMEs),
                                 `Day 1 vs 0`=unique(absTrendDifference[Comparison=="Day 1 vs 0"]),
                                 `Day 2 vs 1`=unique(absTrendDifference[Comparison=="Day 2 vs 1"]),
                                 `Day 3 vs 2`=unique(absTrendDifference[Comparison=="Day 3 vs 2"]),
                                 `Day 4 vs 3`=unique(absTrendDifference[Comparison=="Day 4 vs 3"]),
                                 `Day 5 vs 4`=unique(absTrendDifference[Comparison=="Day 5 vs 4"])),
                              by="Gene.IDs"])

row.names(TABLE_MATRIX)=TABLE_MATRIX$Gene.NAMEs
TABLE_MATRIX=TABLE_MATRIX[,-c(1:2)]

#####################################################################
selection=row.names(TABLE_MATRIX) %in% atLeastOneSigDifference

seed=2
nCluster=5
set.seed(seed)

z_score = function(X){return((X-mean(X))/sd(X))}
M_z=t(apply(TABLE_MATRIX[selection,],1,z_score))

kmeans = pheatmap(M_z,
                  kmeans_k=nCluster,
                  cluster_cols=F,
                  cluster_rows=F
)

bestCluster=apply(kmeans$kmeans$centers,2,which.max)

help=function(X,bestCluster=bestCluster){
  return(which(X==bestCluster))
}
if(length(unique(bestCluster))==length(bestCluster)){
  kmeans_clustering=unlist(lapply(kmeans$kmeans$cluster, help, bestCluster=unname(bestCluster)))
}else{
  kmeans_clustering=kmeans$kmeans$cluster
  
}

M_z_cluster=cbind(M_z,10*kmeans_clustering)

cluster_rows= hclust(dist(M_z_cluster, method="euclidean"),method="average")

annotation_row = data.frame(
  Cluster = as.factor(kmeans_clustering)
)
row.names(annotation_row)=row.names(M_z)

pheatmap(M_z,
         scale = "none",
         cluster_rows=cluster_rows,
         annotation_row = annotation_row,
         cluster_cols=F,
         show_rownames = T,
         color=rev(hcl.colors(10, "Purples 3")),
         cutree_rows = nCluster,
         treeheight_row = 0, 
         filename = paste(prefix,"heatmap.pdf",sep=""))

###############################################################################


