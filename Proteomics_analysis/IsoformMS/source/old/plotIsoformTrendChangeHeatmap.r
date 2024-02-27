list.of.packages <- c("data.table","pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(pheatmap))

MS.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3/de-novo/IsoformTrendchanges.txt"

out="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3/de-novo/pdf/test.pdf"

widthInch=12
fontsize_row=12
fontsize_col=12
cellwidth=10
method="mean"
dayThr=3
cutree_rows=6
distMethod="euclidean"
clusteringMethod="average"
args <- commandArgs(trailingOnly = TRUE) 

MS.path = args[1]
out = args[2]
out2 = args[3]
prefix=args[4]

TABLE = as.data.table(read.table(MS.path, header = T, sep="\t"))

#TODO consider more than two transcripts
TABLE=TABLE[,.(trendDifferences_mean_Day1vsDay0=abs(mean_Day1vsDay0[1]-mean_Day1vsDay0[2]),
               trendDifferences_mean_Day2vsDay1=abs(mean_Day2vsDay1[1]-mean_Day2vsDay1[2]),
               trendDifferences_mean_Day3vsDay2=abs(mean_Day3vsDay2[1]-mean_Day3vsDay2[2]),
               trendDifferences_mean_Day4vsDay3=abs(mean_Day4vsDay3[1]-mean_Day4vsDay3[2]),
               trendDifferences_mean_Day5vsDay4=abs(mean_Day5vsDay4[1]-mean_Day5vsDay4[2]),
               trendDifferences_sum_Day1vsDay0=abs(sum_Day1vsDay0[1]-sum_Day1vsDay0[2]),
               trendDifferences_sum_Day2vsDay1=abs(sum_Day2vsDay1[1]-sum_Day2vsDay1[2]),
               trendDifferences_sum_Day3vsDay2=abs(sum_Day3vsDay2[1]-sum_Day3vsDay2[2]),
               trendDifferences_sum_Day4vsDay3=abs(sum_Day4vsDay3[1]-sum_Day4vsDay3[2]),
               trendDifferences_sum_Day5vsDay4=abs(sum_Day5vsDay4[1]-sum_Day5vsDay4[2])
               ), 
               by=c("gene_id", "gene_name")]

M1=as.matrix(TABLE[,c("trendDifferences_mean_Day1vsDay0",
                      "trendDifferences_mean_Day2vsDay1",
                      "trendDifferences_mean_Day3vsDay2",
                      "trendDifferences_mean_Day4vsDay3",
                      "trendDifferences_mean_Day5vsDay4")])
colnames(M1)=c("Day1\nvs\nDay0","Day2\nvs\nDay1","Day3\nvs\nDay2","Day4\nvs\nDay3","Day5\nvs\nDay4")
row.names(M1)=do.call(rbind,strsplit(TABLE$gene_name,split="\\|"))[,1]
M1=M1[rowSums(is.na(M1))<=dayThr,]

replaceMissingValuesByZero = function(X){
  X[is.na(X)]=0
  return(X)
}

seed=2
removeCluster=2
nCluster=cutree_rows

M1_0=replaceMissingValuesByZero(M1)

set.seed(seed)
kmeans1= pheatmap(M1_0,
                  kmeans_k=nCluster,
                  cluster_cols=F,
                  cluster_rows=F
)

kmeans_clustering1=kmeans1$kmeans$cluster

M1_0=cbind(M1_0,10*kmeans_clustering1)

M2=M1[kmeans_clustering1!=removeCluster,]
M2_0=M1_0[kmeans_clustering1!=removeCluster,]


cluster_rows1= hclust(dist(M1_0, method=distMethod),method=clusteringMethod)
cluster_rows2= hclust(dist(M2_0, method=distMethod),method=clusteringMethod)

clusteN1=kmeans_clustering1
clusteN2=kmeans_clustering1[kmeans_clustering1!=removeCluster]

for (i in 1:nCluster){
  clusteN1[clusteN1==i]=paste(i," (n=",sum(clusteN1==i),")",sep="")  
}

temp=1:(nCluster)
for (i in temp[-c(removeCluster)]){
  clusteN2[clusteN2==i]=paste(i," (n=",sum(clusteN2==i),")",sep="")  
}


annotation_row = data.frame(
  Cluster = as.factor(clusteN1)
)
row.names(annotation_row)=row.names(M1_0)

annotation_row2 = data.frame(
  Cluster = as.factor(clusteN2)
)
row.names(annotation_row2)=row.names(M2_0)

annotation_col = data.frame(
  Day = as.factor(c("day 1 vs. day 0",
                    "day 2 vs. day 1",
                    "day 3 vs. day 2",
                    "day 4 vs. day 3",
                    "day 5 vs. day 4"))
)
row.names(annotation_col)=colnames(M1)

pheatmap(M1,
         cluster_cols=F,
         show_rownames = F,
         cluster_rows= cluster_rows1,
         annotation_row=annotation_row,
         color=rev(hcl.colors(10, "Purples 3")),
         cutree_rows = cutree_rows,
         treeheight_row = 0,
         width = widthInch,
         filename = out
         )


pheatmap(M2,
         show_rownames = T,
         show_colnames = F,
         cluster_cols=F,
         annotation_row=annotation_row2,
         annotation_col=annotation_col,
         cluster_rows= cluster_rows2,
         color=rev(hcl.colors(10, "Purples 3")),
         cutree_rows = cutree_rows-1,
        treeheight_row = 0,
        width = widthInch,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        cellwidth = cellwidth,
        filename=out2
)

for(i in 1:nCluster){
  write.table(row.names(M1)[kmeans_clustering1==1], row.names = F, col.names = F, quote = F, sep = "\t", file = paste(prefix,"cluster_",i,".txt",sep=""))
}

write.table(row.names(M1), row.names = F, col.names = F, quote = F, sep = "\t", file = paste(prefix,"reference.txt",sep=""))