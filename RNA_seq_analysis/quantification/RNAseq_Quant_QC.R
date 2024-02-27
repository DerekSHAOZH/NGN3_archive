library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(matrixStats)
library(GGally)


plot_out_dir <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Plots/' 
samples <- c('day0_1', 'day0_2', 'day0_3', 'day0_4', 'day0_5', 
             'day1_1', 'day1_2', 'day1_3', 'day1_4', 'day1_5',
             'day2_1', 'day2_2', 'day2_3', 'day2_4', 'day2_5',
             'day3_1', 'day3_2', 'day3_3', 'day3_5', 
             'day5_1', 'day5_2', 'day5_4', 'day5_5', 'day5_6')
cond <- c('day0', 'day0', 'day0', 'day0', 'day0', 
          'day1', 'day1', 'day1', 'day1', 'day1',
          'day2', 'day2', 'day2', 'day2', 'day2',
          'day3', 'day3', 'day3', 'day3', 
          'day5', 'day5', 'day5', 'day5', 'day5')

pval <- 0.05
l2fc <- log2(1.5)

####Salmon: reference genome ####

dire<-'/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToRefTrans/'
res_out_dir <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToRefTrans/' 

###gene-level ####
fns <- paste0(dire, samples,'/quant.genes.sf')


txi<-tximport(fns, type="salmon", txIn = FALSE, txOut = FALSE, geneIdCol = 'Name')
txi$length[txi$length == 0] <- 1
colnames(txi$counts) <- samples

#export gene-level TPM abundance
# gene_TPM <- txi$abundance
# colnames(gene_TPM) <- samples
# write.table(gene_TPM, paste0(res_out_dir, "salmon_gene_TPM.txt"), quote = F, sep = '\t', row.names = T)


df <- data.frame(condition = factor(cond))
rownames(df) <- samples
dds <- DESeqDataSetFromTximport(txi, df, ~condition)
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds))
dds <-  estimateDispersions(dds)
dds <-  nbinomWaldTest(dds)
vsd <- vst(dds)



#heathmap of correlation all samples 
sampleDists <- dist(t(assay(vsd)))
sampleDists <- cor(assay(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$condition
corrMatrix<-sampleDistMatrix
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pdf(paste0(plot_out_dir, 'salmon_correlation_heatmap.pdf'))
p<-pheatmap(corrMatrix,col=colors)
print(p)
dev.off()

#PCA all samples 
pdf(paste0(plot_out_dir,'salmon_PCA.pdf'))
plotPCA(vsd)
dev.off()

#QC
df_counts <- counts(dds)
colnames(df_counts)<- samples
pdf(paste0(plot_out_dir,'salmon_counts_distributions.pdf'))
p<-boxplot(log10(counts(dds)+1))
print(p)
dev.off()





###isoform-level ####
fns <- paste0(dire, samples,'/quant.sf')

txi<-tximport(fns, type="salmon", txIn = TRUE, txOut = TRUE)
txi$length[txi$length == 0] <- 1

transcript_TPM <- txi$abundance
colnames(transcript_TPM) <- samples
write.table(transcript_TPM, paste0(res_out_dir, "salmon_transcript_TPM.txt"), quote = F, sep = '\t', row.names = T)

###compare to RSEM results ####
###gene-level ####
##salmon
dire<-'/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToRefTrans/'
fns <- paste0(dire, samples,'/quant.genes.sf')


txi<-tximport(fns, type="salmon", txIn = FALSE, txOut = FALSE, geneIdCol = 'Name')
txi$length[txi$length == 0] <- 1

salmon_TPM <- txi$abundance
colnames(salmon_TPM) <- samples

#filter by median TPM > 1
salmon_TPM <- salmon_TPM[rowMedians(salmon_TPM) > 1, ]
#filter out "PAR_Y" gene
salmon_TPM <- salmon_TPM[!endsWith(rownames(salmon_TPM), "PAR_Y"), ]

##rsem
rsem_dire<-'/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/'

rsem_samples <- c('OJ63','OJ64','OJ65','OJ66','OJ67','OJ68','OJ69','OJ70','OJ71','OJ72','OJ73','OJ74',
             'OJ75','OJ76','OJ77','OJ79','OJ80','OJ81','OJ83','OJ85','OJ86','OJ87','OJ88','OJ89')

fns <- paste0(rsem_dire,rsem_samples,'.genes.results')

txi<-tximport(fns, type="rsem", txIn = FALSE, txOut = FALSE)
txi$length[txi$length == 0] <- 1

rsem_TPM <- txi$abundance
colnames(rsem_TPM) <- samples

#filter by median TPM > 1
rsem_TPM <- rsem_TPM[rowMedians(rsem_TPM) > 1, ]
#filter out "PAR_Y" gene
rsem_TPM <- rsem_TPM[!endsWith(rownames(rsem_TPM), "PAR_Y"), ]

##common gene: 5197 genes
common <- intersect(rownames(salmon_TPM), rownames(rsem_TPM))

salmon_TPM <- salmon_TPM[common, ]
rsem_TPM <- rsem_TPM[common, ]

corr <- cor(salmon_TPM, rsem_TPM, method = "spearman")
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pdf(paste0(plot_out_dir,'salmon_rsem_gene_TPM_correlation_heatmap.pdf'))
pheatmap(corr,col=colors
         , cluster_rows = F, cluster_cols = F,
)
dev.off()

##Pairwise comparisons

#Spearman's correlations
gene_corr <- diag(corr)
#max: 0.9726608
max(gene_corr)
#min: 0.9471051
min(gene_corr)
#median: 0.9592002
median(gene_corr)


#Scatter plot 
#log2 transformation with a pseudo-count of 1
salmon_log2_TPM <- log2(salmon_TPM + 1)
rsem_log2_TPM <- log2(rsem_TPM + 1)


pdf(paste0(plot_out_dir,'salmon_rsem_gene_TPM_correlation_scatterplot.pdf'))
for (x in samples){
  df <- data.frame("Salmon" = salmon_log2_TPM[, x],
                   "RSEM" = rsem_log2_TPM[, x])
  
  p <- ggpairs(df, upper = list(continuous = wrap("cor", method = "spearman", size = 8)),
          title = x) + 
    theme(axis.text = element_text(size = 20),
          strip.text = element_text(size = 20),
          plot.title = element_text(size=20))
  print(p)
}

dev.off() 

###isoform-level ####
##salmon
fns <- paste0(dire, samples,'/quant.sf')


txi<-tximport(fns, type="salmon", txIn = TRUE, txOut = TRUE)
txi$length[txi$length == 0] <- 1

salmon_TPM <- txi$abundance
colnames(salmon_TPM) <- samples

#filter by median TPM > 1
salmon_TPM <- salmon_TPM[rowMedians(salmon_TPM) > 1, ]
#filter out "PAR_Y" isoforms
salmon_TPM <- salmon_TPM[!endsWith(rownames(salmon_TPM), "PAR_Y"), ]

##rsem
rsem_dire<-'/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/'

rsem_samples <- c('OJ63','OJ64','OJ65','OJ66','OJ67','OJ68','OJ69','OJ70','OJ71','OJ72','OJ73','OJ74',
                  'OJ75','OJ76','OJ77','OJ79','OJ80','OJ81','OJ83','OJ85','OJ86','OJ87','OJ88','OJ89')

fns <- paste0(rsem_dire,rsem_samples,'.isoforms.results')

txi<-tximport(fns, type="rsem", txIn = TRUE, txOut = TRUE)
txi$length[txi$length == 0] <- 1

rsem_TPM <- txi$abundance
colnames(rsem_TPM) <- samples

#filter by median TPM > 1
rsem_TPM <- rsem_TPM[rowMedians(rsem_TPM) > 1, ]
#filter out "PAR_Y" isoforms
rsem_TPM <- rsem_TPM[!endsWith(rownames(rsem_TPM), "PAR_Y"), ]

##common isoforms: 31104 isoforms
common <- intersect(rownames(salmon_TPM), rownames(rsem_TPM))

salmon_TPM <- salmon_TPM[common, ]
rsem_TPM <- rsem_TPM[common, ]

corr <- cor(salmon_TPM, rsem_TPM, method = "spearman")
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pdf(paste0(plot_out_dir,'salmon_rsem_isoform_TPM_correlation_heatmap.pdf'))
pheatmap(corr,col=colors
         , cluster_rows = F, cluster_cols = F,
)
dev.off()

##Pairwise comparisons

#Spearman's correlations
isoform_corr <- diag(corr)
#max: 0.9120054
max(isoform_corr)
#min: 0.8879579
min(isoform_corr)
#median: 0.8979895
median(isoform_corr)


#Scatter plot 
#log2 transformation with a pseudo-count of 1
salmon_log2_TPM <- log2(salmon_TPM + 1)
rsem_log2_TPM <- log2(rsem_TPM + 1)


pdf(paste0(plot_out_dir,'salmon_rsem_isoform_TPM_correlation_scatterplot.pdf'))
for (x in samples){
  df <- data.frame("Salmon" = salmon_log2_TPM[, x],
                   "RSEM" = rsem_log2_TPM[, x])
  
  p <- ggpairs(df, upper = list(continuous = wrap("cor", method = "spearman", size = 8)),
               title = x) + 
    theme(axis.text = element_text(size = 20),
          strip.text = element_text(size = 20),
          plot.title = element_text(size=20))
  print(p)
}

dev.off() 

###boxplot of correlations ####
df <- data.frame(corr = c(gene_corr, isoform_corr), group = c(rep("Genes\nn=5197", 24), rep("Transcript isoforms\nn=31104", 24)))

pdf(file = paste0(plot_out_dir, "salmon_rsem_TPM_correlation_boxplot.pdf"), 4,4)
ggplot(df, aes(x = group, y = corr, fill = group))+
  geom_boxplot(alpha = 0.8)+
  # geom_violin(alpha = 0.5)+ 
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              dotsize = 0.5) +
  
  geom_jitter(position = position_jitter(seed = 1, width = 0.3)) +
  ylab("Spearman's correlation coefficient") + 
  xlab("Targets of expression quantification") + 
  
  
  theme_classic() +
  theme(legend.position="none", text = element_text(size = 15))
dev.off()



####RSEM: ONT transcriptome ####
dire<-'/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/RSEM_Quantification/'
res_out_dir <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/RSEM_Quantification/' 
DTE_res_out_dir <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/DTE_analysis_with_RSEM/' 

###isoform-level ####
fns <- paste0(dire, samples,'/', samples, '.isoforms.results')

txi<-tximport(fns, type="rsem", txIn = TRUE, txOut = TRUE)
txi$length[txi$length == 0] <- 1

# transcript_TPM <- txi$abundance
# colnames(transcript_TPM) <- samples
# write.table(transcript_TPM, paste0(res_out_dir, "rsem_transcript_TPM.txt"), quote = F, sep = '\t', row.names = T)

df <- data.frame(condition = factor(cond))
rownames(df) <- colnames(txi$counts)


dds <- DESeqDataSetFromTximport(txi, df, ~condition)

#QC ####
#sample correlation
vsd <- vst(dds)
sampleDists <- dist(t(assay(vsd)))
sampleDists <- cor(assay(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- cond
colnames(sampleDistMatrix) <- cond
corrMatrix<-sampleDistMatrix
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pdf(paste0(plot_out_dir,'rsem_OntTrans_correlation_heatmap.pdf'))
pdf(paste0(plot_out_dir,'rsem_OntTrans_correlation_heatmap_unclustered.pdf'))
p<-pheatmap(corrMatrix,col=colors
            , cluster_rows = F, cluster_cols = F,
            )
print(p)
dev.off()


#count distribution 
df_counts <- counts(dds)
colnames(df_counts)<- samples
pdf(paste0(plot_out_dir,'rsem_OntTrans_counts_distributions.pdf'))
boxplot(log10(counts(dds)+1),names=samples)
dev.off()

#PCA of samples
vsd <- vst(dds)
pdf(paste0(plot_out_dir,'rsem_OntTrans_transcript_PCA.pdf'))
plotPCA(vsd)
dev.off()

#finding DE transcript ####
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds))

Pcounts <- counts(dds, normalized=TRUE)
colnames(Pcounts) <- samples
write.csv(Pcounts, '/project/Neurodifferentiation_System/owlmayerTemporary/derek/RNAseq_quantification/Results/AlignedToOntTrans/DTE_analysis_with_RSEM/count_transcript_deseq2norm.csv', quote = F)

dds <-  estimateDispersions(dds)
dds <-  nbinomWaldTest(dds)



deset <-c()
for (i in 1:(length(levels(df$condition))-1)){
  for (j in (i+1):length(levels(df$condition))){
    control<-levels(df$condition)[i]
    condition<-levels(df$condition)[j]
    dir.create(paste0(DTE_res_out_dir,condition,'_background_',control,'/'))
    
    res <- results(dds,contrast=c("condition",condition,control),alpha=pval,lfcThreshold=l2fc)
    res <- lfcShrink(dds, contrast=c("condition",condition,control), res=res,type="ashr")
    
    results <- as.data.frame(res)
    deid <- results %>% filter(padj<pval) %>% filter(abs(log2FoldChange)>l2fc) %>% rownames()
    deset<- c(deset,deid)
    print(paste(control,condition,length(deid),length(deset)))
    results$transcript_id<-rownames(res)
    # results<-results %>% separate(ensgene,into='ensgene',sep='\\.')
    # results<-results%>% left_join(grch38[, c("ensgene", "symbol","biotype", "description")], by = "ensgene")
    results %>% write.table(paste0(DTE_res_out_dir,condition,'_background_',control,'/DEseq_all_res.csv'), quote=F, row.names = F)
    
    #creating tables for GO
    go_lists <- results %>% filter(baseMean>0)
    rownames(go_lists) <- c()
    
    go_lists %>% select(transcript_id) %>%
      write.table(paste0(DTE_res_out_dir,condition,'_background_',control,'/background_list.txt'),row.names=F,col.names=F,quote=F)
    
    go_lists %>%
      filter(padj<pval) %>%
      filter(log2FoldChange>l2fc) %>%
      select(transcript_id) %>%
      write.table(paste0(DTE_res_out_dir,condition,'_background_',control,'/upregulated_list.txt'),row.names=F,col.names=F,quote=F)
    
    go_lists %>%
      filter(padj<pval) %>%
      filter(log2FoldChange<l2fc) %>%
      select(transcript_id) %>%
      write.table(paste0(DTE_res_out_dir,condition,'_background_',control,'/downregulated_list.txt'),row.names=F,col.names=F,quote=F)
    
    #saving differentially expressed genes
    sigres <- results %>%
      filter(padj<pval) %>%
      filter(baseMean>0)%>%
      filter(abs(log2FoldChange)>l2fc)%>%
      arrange(desc(abs(log2FoldChange)),desc(baseMean))
    sigres %>% write.csv(file=paste0(DTE_res_out_dir,condition,'_background_',control,'/diff_expression.csv'),quote=F, row.names = F)
    
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
    pdf(paste0(DTE_res_out_dir,condition,'_background_',control,'/Volcano_plot.pdf'))
    plot(results$log2FoldChange,-1*log10(results$padj),pch = 20,
         xlab = "Log2FoldChange", ylab = "-log10(padj)", col = col.vec, xlim = c(-1*x_lim, x_lim))
    legend("topright", paste(c("up = ","down = "),c(n.differential.up, n.differential.down), sep = ""), pch = 20, col = c(red,blue))
    dev.off()
  }
}

de_genes <- unique(deset)
write(de_genes, paste0(DTE_res_out_dir,'DE_at_any_timepoint.txt'))
