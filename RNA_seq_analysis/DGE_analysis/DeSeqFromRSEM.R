#loading the packages
library(tidyverse)
library(tximport)
library(DESeq2)
library(annotables)
library(ape)
library(scico)
library(colorspace)
library(RColorBrewer)
library(pheatmap)

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

#loading the RSEM results, creating the metadata table
txi<-tximport(fns, type="rsem", txIn = FALSE, txOut = FALSE)
txi$length[txi$length == 0] <- 1
df <- data.frame(condition = factor(cond))
rownames(df) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, df, ~condition)

vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDists <- cor(assay(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- cond
colnames(sampleDistMatrix) <- cond
corrMatrix<-sampleDistMatrix
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pdf(paste0(direout,'sample_correlations.pdf'))
p<-pheatmap(corrMatrix,col=colors)
print(p)
dev.off()

#QC
df_counts <- counts(dds)
colnames(df_counts)<- samples
pdf(paste0(direout,'counts_distributions.pdf'))
boxplot(log10(counts(dds)+1),names=samples)
dev.off()

#PCA of samples
vsd <- vst(dds)
pdf(paste0(direout,'sample_PCA.pdf'))
plotPCA(vsd)
dev.off()

#finding DE genes
View(counts(dds))
dds <- estimateSizeFactors(dds)
normalizationFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)




colnames(normalized_counts) <- paste0(cond, '_', sapply(gsub("rep","", rep), function(x) (x)))
View(normalized_counts)
write.table(normalized_counts, file=paste0(direout, "counts_normalized.txt"), sep="\t", quote=F, col.names=)




sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds))
dds <-  estimateDispersions(dds)
dds <-  nbinomWaldTest(dds)

deset <-c()
for (i in 1:(length(levels(df$condition))-1)){
  for (j in (i+1):length(levels(df$condition))){
    control<-levels(df$condition)[i]
    condition<-levels(df$condition)[j]
    dir.create(paste0(direout,condition,'_background_',control,'/'))
    
    res <- results(dds,contrast=c("condition",condition,control),alpha=pval,lfcThreshold=l2fc)
    res <- lfcShrink(dds, contrast=c("condition",condition,control), res=res,type="ashr")
    
    results <- as.data.frame(res)
    deid <- results %>% filter(padj<pval) %>% filter(abs(log2FoldChange)>l2fc) %>% rownames()
    deset<- c(deset,deid)
    print(paste(control,condition,length(deid),length(deset)))
    results$ensgene<-rownames(res)
    results<-results %>% separate(ensgene,into='ensgene',sep='\\.')
    results<-results%>% left_join(grch38[, c("ensgene", "symbol","biotype", "description")], by = "ensgene")
    results %>% write.table(paste0(direout,condition,'_background_',control,'/DEseq_all_res.csv'), quote=F)
    
    #creating tables for GO
    go_lists <- results %>% filter(baseMean>0)
    rownames(go_lists) <- c()
    
    go_lists %>% select(ensgene) %>%
      write.table(paste0(direout,condition,'_background_',control,'/background_list.txt'),row.names=F,col.names=F,quote=F)
    
    go_lists %>%
      filter(padj<pval) %>%
      filter(log2FoldChange>l2fc) %>%
      select(ensgene) %>%
      write.table(paste0(direout,condition,'_background_',control,'/upregulated_list.txt'),row.names=F,col.names=F,quote=F)
    
    go_lists %>%
      filter(padj<pval) %>%
      filter(log2FoldChange<l2fc) %>%
      select(ensgene) %>%
      write.table(paste0(direout,condition,'_background_',control,'/downregulated_list.txt'),row.names=F,col.names=F,quote=F)
    
    #saving differentially expressed genes
    sigres <- results %>%
      filter(padj<pval) %>%
      filter(baseMean>0)%>%
      filter(abs(log2FoldChange)>l2fc)%>%
      arrange(desc(abs(log2FoldChange)),desc(baseMean))
      sigres %>% write.csv(file=paste0(direout,condition,'_background_',control,'/diff_expression.csv'))
    
    #Numbers of differentially expressed genes
    n.differential = nrow(sigres)
    n.differential.up = sum(sigres$log2FoldChange > 0)
    n.differential.down = sum(sigres$log2FoldChange < 0)
    print(n.differential)
    print(n.differential.up)
    print(n.differential.down)
    
    #Volcano plot
    grey = rgb(51, 51, 51,  50, maxColorValue = 255)
    red = rgb(255, 0, 0, 50, maxColorValue = 255)
    blue = rgb(0, 0, 255, 50, maxColorValue = 255)
    col.vec = rep(grey,nrow(res))
    col.vec[results$padj < pval & results$log2FoldChange < l2fc] = blue
    col.vec[results$padj < pval & results$log2FoldChange > l2fc] = red
    x_lim = max(abs(results$log2FoldChange), na.rm = T)
    pdf(paste0(direout,condition,'_background_',control,'/Volcano_plot.pdf'))
    plot(results$log2FoldChange,-1*log10(results$padj),pch = 20,
         xlab = "Log2FoldChange", ylab = "-log10(padj)", col = col.vec, xlim = c(-1*x_lim, x_lim))
    legend("topright", paste(c("up = ","down = "),c(n.differential.up, n.differential.down), sep = ""), pch = 20, col = c(red,blue))
    dev.off()
  }
}

# de_genes <- unique(deset)
# normCounts<-counts(dds,normalized=T)
# DEnormCounts<-normCounts[de_genes,]
# colnames(DEnormCounts)<-paste(cond,rep,sep='_')
# dfDEcounts<- as.data.frame(DEnormCounts)
# 
# for (c in cond){
#   dfDEcounts[,paste0("Mean_",c)]<-dfDEcounts%>% select(starts_with(c))%>%rowMeans()
# }
# meanCounts<-dfDEcounts%>% select(starts_with("Mean"))
# standMeanCounts<-as.data.frame(t(scale(t(meanCounts))))
# kclusters<-kmeans(standMeanCounts,7)
# standMeanCounts$class<-kclusters$cluster
# standMeanCounts<- standMeanCounts %>% arrange(class)
# kclasses<-standMeanCounts$class
# standMeanCounts$class<-NA
# pheatmap(standMeanCounts,show_rownames=F,cluster_cols=F,cluster_rows = F)

de_genes <- unique(deset)
write(de_genes, paste0(direout,'DE_at_any_timepoint.txt'))
deTPM<-txi$abundance[de_genes,]
colnames(deTPM)<-paste(cond,rep,sep='_')
deTPM<- as.data.frame(deTPM)

for (c in cond){
  deTPM[,paste0("Mean_",c)]<-deTPM%>% select(starts_with(c))%>%rowMeans()
}
meanTPM<-deTPM%>% select(starts_with("Mean"))
standTPM<-as.data.frame(t(scale(t(meanTPM))))
kclusters<-kmeans(standTPM,4)
standTPM$max<- sapply(as.data.frame(t(standTPM)),max)
standTPM$maxi<-max.col(standTPM)
standTPM$class<-kclusters$cluster
standTPM<- standTPM %>% arrange(class)
kclasses<-standTPM %>%select(class)
standTPM<-standTPM%>%select(-c(class,max,maxi))

pdf(paste0(direout,'Heatmap.pdf'))
pheatmap(standTPM,show_rownames=F,cluster_cols=F,cluster_rows = F, annotation_colors=list("class"=qualitative_hcl(4,palette="Dark3")),
         annotation_row=kclasses,color =rev(scico(100,palette='roma'))[10:90])
dev.off()


meanTPM$ensgene <- row.names(meanTPM)
meanTPM <- meanTPM %>% separate(ensgene,c('ensgene'))
meanTPM <- meanTPM %>% left_join(grch38,by='ensgene')
meanTPM %>% write.table(paste0(direout,'/DEgenes_info.csv'), quote=F,sep='\t',row.names = F)