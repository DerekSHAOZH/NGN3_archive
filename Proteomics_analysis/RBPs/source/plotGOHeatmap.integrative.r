list.of.packages <- c("org.Hs.eg.db","GO.db","pheatmap","data.table","scico")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(scico))
suppressPackageStartupMessages(library(data.table))


GO_term="GO:0008380;GO:0003723"
TABLE.ONT.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/ONT-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt"
TABLE.RNA.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/RNA-seq_NeuroDifferentation_NGN3/gene/countTable.RLE.txt"

significant.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/RBPs/significant.txt"
nCluster=3
out_reference="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/RBPs/RBP_splicing.txt"
out="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/RBPs/RBP_splicing.pdf"

args <- commandArgs(trailingOnly = TRUE) 



GO_term=args[1]
TABLE.ONT.path=args[2]
TABLE.RNA.path=args[3]
significant.path=args[4]
nCluster=as.numeric(args[5])
out_GO=args[6]
out=args[7]
out2=args[8]
out_cluster=args[9]
out_reference=args[10]


GO=unlist(strsplit(GO_term,split=";"))
significant =read.table(significant.path, header=F)$V1

ANNOTATION_TABLE=AnnotationDbi::select(org.Hs.eg.db, keys=GO, columns = c("ENSEMBL","SYMBOL"), keytype = "GOALL")


removeManually=function(X){
  # remove wrong entries manually
  return(X[!(
    (X$ENSEMBL=="ENSG00000205571" & X$SYMBOL=="SMN1") |
    (X$ENSEMBL=="ENSG00000266086" & X$SYMBOL=="SRSF1") |
    (X$ENSEMBL=="ENSG00000280987" & X$SYMBOL=="MATR3") |
    (X$ENSEMBL=="ENSG00000276345" & X$SYMBOL=="MRPL23") |
    (X$ENSEMBL=="ENSG00000288141" & X$SYMBOL=="MRPL23") |
    (X$ENSEMBL=="ENSG00000284764" & X$SYMBOL=="PINX1") |
    (X$ENSEMBL=="ENSG00000258724" & X$SYMBOL=="PINX1") |
      (X$ENSEMBL=="ENSG00000275740" & X$SYMBOL=="RBM27") |
      (X$ENSEMBL=="ENSG00000255508" & X$SYMBOL=="TUT1") #|
      #(X$ENSEMBL=="ENSG00000284764" & X$SYMBOL=="PINX1") |
      #(X$ENSEMBL=="ENSG00000284764" & X$SYMBOL=="PINX1") |
  ),])
}

RBPs=removeManually(na.omit(unique(ANNOTATION_TABLE[,c("ENSEMBL","SYMBOL")])))

write.table(unique(RBPs), row.names = F, col.names = c("ENSEMBL ID","Gene name"), quote=F, sep="\t", file=out_GO)

removeManually2=function(X){
  return(X[!(X$ENSEMBL=="ENSG00000273784" & X$SYMBOL=="HNRNPA1L2"),])
}

TABLE_ONT=read.table(TABLE.ONT.path, header=T)
TABLE_RNA=read.table(TABLE.RNA.path, header=T)
TABLE_RNA=TABLE_RNA[,grepl("day[0,3,5]",colnames(TABLE_RNA))]

TABLE_ONT$ENSEMBL=unlist(do.call(rbind,strsplit(row.names(TABLE_ONT),split="\\."))[,1])
TABLE_RNA$ENSEMBL=unlist(do.call(rbind,strsplit(row.names(TABLE_RNA),split="\\."))[,1])

TABLE=merge(TABLE_ONT, TABLE_RNA, by="ENSEMBL",all=T)

TABLE=unique(merge(TABLE,RBPs,by="ENSEMBL",all.x=T))
TABLE=removeManually2(TABLE)


TABLE_splicing_factors=TABLE[TABLE$ENSEMBL%in%RBPs$ENSEMBL,]
row.names(TABLE_splicing_factors)=TABLE_splicing_factors$SYMBOL

TABLE_splicing_factors_significant=TABLE[TABLE$ENSEMBL%in%RBPs$ENSEMBL & TABLE$ENSEMBL%in%significant,]
row.names(TABLE_splicing_factors_significant)=TABLE_splicing_factors_significant$SYMBOL

M=TABLE_splicing_factors_significant[,grepl("day",colnames(TABLE_splicing_factors_significant))]


# replace NA as 0
M_0=M
M_0[is.na(M_0)]=0
z_score = function(X){
  if(sum(X)==0 || sum(is.na(X))!=0){
    return(X)
  }else{
    return((X-mean(X))/sd(X))
  }
}
# 
M_z=cbind(t(apply(M_0[,grepl("ONT_",colnames(M_0))],1,z_score)),
          t(apply(M_0[,grepl("RNA_",colnames(M_0))],1,z_score)))

# visualization matrix
M_visualize=cbind(t(apply(M[,grepl("ONT_",colnames(M))],1,z_score)),
          t(apply(M[,grepl("RNA_",colnames(M))],1,z_score)))

seed=0
set.seed(seed)


kmeans_1= pheatmap(M_z,
                   scale="none",
                   kmeans_k= nCluster,
                   cluster_cols=F,
                   cluster_rows=F
)


kmeans_clustering_1=kmeans_1$kmeans$cluster

TABLE_splicing_factors_significant_cluster=cbind(M_z, 10*kmeans_clustering_1)
cluster_rows_1= hclust(dist(TABLE_splicing_factors_significant_cluster))

###

annotation_row = data.frame(
  Cluster=as.factor(kmeans_clustering_1)
)
row.names(annotation_row)=row.names(TABLE_splicing_factors_significant)

annotation_col = data.frame(
  Replicate = as.factor(do.call(rbind,strsplit(colnames(M),split="_"))[,3]),
  Day = as.factor(substring(as.factor(do.call(rbind,strsplit(colnames(M),split="_"))[,2]),4)),
  Method = as.factor(do.call(rbind,strsplit(colnames(M),split="_"))[,1])

)
row.names(annotation_col)=colnames(M)

pheatmap(M_visualize[,grepl("ONT_",colnames(M))],
         scale="none",
         cluster_cols=F,
         show_rownames = F,
         cluster_rows= cluster_rows_1,
         annotation_col=annotation_col,
         show_colnames = F,
         annotation_row=annotation_row,
         color=rev(scico(100,palette='roma'))[10:90],
         cutree_rows = nCluster,
         treeheight_row = 0 ,
         filename = out
)

pheatmap(M_visualize[,grepl("RNA_",colnames(M))],
         scale="none",
         cluster_cols=F,
         show_rownames = F,
         cluster_rows= cluster_rows_1,
         annotation_col=annotation_col,
         show_colnames = F,
         annotation_row=annotation_row,
         color=rev(scico(100,palette='roma'))[10:90],
         cutree_rows = nCluster,
         treeheight_row = 0,
         filename = out2
)



names(kmeans_clustering_1)=TABLE_splicing_factors_significant$ENSEMBL

write.table(sort(kmeans_clustering_1), col.names = F, quote = F, sep = "\t", file=out_cluster)
write.table(unique(TABLE_splicing_factors$ENSEMBL), row.names=F, col.names = F, quote = F, sep = "\t", file=out_reference)
