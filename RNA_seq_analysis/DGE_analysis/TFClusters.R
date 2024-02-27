#loading the packages
library(tidyverse)
library(tximport)
library(DESeq2)
library(annotables)
library(ape)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(scico)
library(colorspace)
library(RColorBrewer)
library(pheatmap)
library(conflicted)
library(rtracklayer)
library(scico)
library(cluster)
library(factoextra)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

#reading in the paths and metadata
annotfn <- '/project/PausingDynamics/GeneralResources/Genecode29/gencode.v29.annotation.sorted.gff3'
dire<-'/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/'
direout<-'/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/'

samples <- c('OJ63','OJ64','OJ65','OJ66','OJ67','OJ68','OJ69','OJ70','OJ71','OJ72','OJ73','OJ74',
             'OJ75','OJ76','OJ77','OJ79','OJ80','OJ81','OJ83','OJ85','OJ86','OJ87','OJ88','OJ89')
cond <- c('day0','day0','day0','day0','day0','day1','day1','day1','day1','day1','day2','day2',
          'day2','day2','day2','day3','day3','day3','day3','day5','day5','day5','day5','day5')
rep <- c('rep1','rep2','rep3','rep4','rep5','rep1','rep2','rep3','rep4','rep5','rep1','rep2',
          'rep3','rep4','rep5','rep1','rep2','rep3','rep5','rep1','rep2','rep4','rep5','rep6')
fns <- paste0(dire,samples,'.genes.results')
pval <- 0.05
l2fc <- log2(1.5)

qvalgo<-0.05

colors <- data.frame("low" = c("#c5a4a4","#a4c5ab","#a4aac5"), 
                     "high" = c("#ff0000","#6aff89","#6a87ffff"),
                     "symbol"=c("BP","CC","MF"))
rownames(colors) <- c("biological process","cellular compartment","molecular function")

#loading the RSEM results, creating the metadata table
txi<-tximport(fns, type="rsem", txIn = FALSE, txOut = FALSE)
txi$length[txi$length == 0] <- 1
df <- data.frame(condition = factor(cond))
rownames(df) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, df, ~condition)
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds))
dds <-  estimateDispersions(dds)
dds <-  nbinomWaldTest(dds)

deset <-c()
for (i in 1:(length(levels(df$condition))-1)){
  for (j in (i+1):length(levels(df$condition))){
    control<-levels(df$condition)[i]
    condition<-levels(df$condition)[j]
    res <- results(dds,contrast=c("condition",condition,control),alpha=pval,lfcThreshold=l2fc)
    res <- lfcShrink(dds, contrast=c("condition",condition,control), res=res,type="ashr")
    results <- as.data.frame(res)
    deid <- results %>% filter(padj<pval) %>% filter(abs(log2FoldChange)>l2fc) %>% rownames()
    deset<- c(deset,deid)
  }
}

de_genes <- unique(deset)
deTPM<-txi$abundance[de_genes,]
colnames(deTPM)<-paste(cond,rep,sep='_')
deTPM<- as.data.frame(deTPM)

TF_fn<-'/project/Neurodifferentiation_System/GeneralResources/List_of_genes/TF/human_TF.txt'
tfs<-read.csv(TF_fn, col.names = 'TF')
tfs<-as.vector(as.character(map(strsplit(tfs$TF,'\\.'),1)))
deTPM$ID<-as.vector(as.character(map(strsplit(rownames(deTPM),'\\.'),1)))
deTF<-deTPM %>% filter(ID %in% tfs)

for (c in cond){
  deTF[,paste0("Mean_",c)]<-deTF %>% select(starts_with(c))%>%rowMeans()
}
meanTF<-deTF%>% select(starts_with("Mean"))
standTF<-as.data.frame(t(scale(t(meanTF))))

df <- read.csv("/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/TF_clusters.csv", row.names = 1)
df <- df[!duplicated(df$symbol), ]
df$class <- as.numeric(df$class)
df <- df[order(df$class, decreasing = T), ]

rownames(df) <- df$symbol
standTF <- df%>% select(starts_with("Mean"))

# kclasses<-standTF %>%select(class)

# kclusters<-kmeans(standTF,3)
# standTF$class<-kclusters$cluster
# standTF<- standTF %>% arrange(class)


kclasses <- df%>%select(class)
kclasses$class <- factor(kclasses$class, levels = c(3,2,1))


# kclasses<-standTF %>%select(class)
# standTF<-standTF%>%select(-c(class))

Var1        <- qualitative_hcl(3,palette="Dark3")
names(Var1) <- factor(c(3,2,1))

pdf(paste0(direout,'Heatmap_TF_derek.pdf'), 6,4)
pheatmap(standTF,show_rownames=F,cluster_cols=F,cluster_rows = F, annotation_colors=list("class"=Var1),
         annotation_row=kclasses,color =rev(scico(100,palette='roma'))[10:90],
         annotation_names_row = F, labels_col = c('day 0', 'day 1', 'day 2',' day 3', 'day 5'), angle_col = '0',
         annotation_legend = F, fontsize = 12, border_color = NA)
dev.off()


