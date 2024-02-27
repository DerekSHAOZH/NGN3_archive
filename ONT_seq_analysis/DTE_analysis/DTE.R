library(dplyr)
library(tidyr)
library(ggplot2)
library(venn)
library(tximport)
library(patchwork)
library(DESeq2)
library(colorspace)
library(RColorBrewer)
library(pheatmap)
library(rtracklayer)
library(scico)
plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Plots/'
res_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Results/DTE_analysis/'
# annotfn <- '/project/Neurodifferentiation_System/GeneralResources/gencode.v29.annotation.sorted.gff3'
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pval <- 0.05
l2fc <- log2(1.5)

###Data preparation ####
#load gene annotation file ####


# annot <- import.gff3(annotfn)
# annot=as.data.frame(annot)
# annot <- annot %>% filter(type == 'gene')
# annot <- annot %>% select(gene_id, gene_name, gene_type)
# annot <- annot[-which(duplicated(annot$gene_id)), ]
# annot <- annot %>% separate(gene_id, into='gene_id', sep='\\.')

#load ONT transcript count data ####
samples <- c('day0_1', 'day0_2', 'day0_3', 'day3_1', 'day3_2', 'day3_3', 'day5_1', 'day5_2', 'day5_3')

#old quantification 
# data_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts.txt'

# data <- read.csv(file=data_fn, header=TRUE, sep=",")
# data <- data[colnames(data)[c(2:(length(samples) + 1), which(colnames(data) == 'transcript_id') )]]
# 
# data <- data[, sort(colnames(data))]


#updated quantification with filtered annotation
data_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Results/counts.txt'

data <- read.table(file=data_fn, header=TRUE, sep="\t")
# data <- data[colnames(data)[c(2:(length(samples) + 1), which(colnames(data) == 'transcript_id') )]]

# data <- data[, sort(colnames(data))]

#remove duplicated transcript ids
data <- data[!duplicated(data$transcript_id), ]


# countData <- aggregate(. ~ gene_id, data = data, FUN = sum)
# countData <- countData[-1, ]
rownames(data) <- data$transcript_id

#for old
# data <- data[, -ncol(data)]
#for updated
data <- data[, -1]

# data <- data[-which(grepl("PAR_Y", rownames(countData), fixed = T)), ]


countData <- as.matrix(data)



df <- data.frame(name = colnames(countData))
df$condition = factor(sapply(strsplit(df$name,"_"), `[`, 1))

###Perform DE analysis ####
dds <- DESeqDataSetFromMatrix(countData,df,design=~condition)
dds <- estimateSizeFactors(dds)
Pcounts <- counts(dds, normalized=TRUE)
# write.csv(Pcounts, '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/DE_analysis/count_deseq2norm.csv', quote = F)
write.csv(Pcounts, '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Results/DTE_analysis/count_transcript_deseq2norm.csv', quote = F)

vsd <- vst(dds)



#heathmap of correlation all samples 
sampleDists <- dist(t(assay(vsd)))
sampleDists <- cor(assay(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$condition
corrMatrix<-sampleDistMatrix
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
pdf(paste0(plot_dir_out, 'transcript_correlation_heatmap.filtered_quantification.pdf'))
pheatmap(corrMatrix,col=colors)
dev.off()

#PCA all samples 
pdf(paste0(plot_dir_out,'transcript_PCA.filtered_quantification.pdf'))
plotPCA(vsd)
dev.off()

#QC
df_counts <- counts(dds)
colnames(df_counts)<- samples
pdf(paste0(plot_dir_out,'transcript_counts_distributions.filtered_quantification.pdf'))
p<-boxplot(log10(counts(dds)+1))
print(p)
dev.off()







#find DE genes ####
sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds))
dds <-  estimateDispersions(dds)
dds <-  nbinomWaldTest(dds)


deset <- c()
for (i in 1:(length(levels(df$condition))-1)){
  for (j in (i+1):length(levels(df$condition))){ 
    control<-levels(df$condition)[i]
    condition<-levels(df$condition)[j]
    dir.create(paste0(res_dir_out,condition,'_vs_',control,'/'))
    
    res <- results(dds,contrast=c("condition",condition,control),alpha=pval,lfcThreshold=l2fc)
    res <- lfcShrink(dds, contrast=c("condition",condition,control), res=res,type="ashr")
    
    results <- as.data.frame(res)
    
    results$transcript_id<-rownames(res)
    # results<-results %>% separate(gene_id,into='gene_id',sep='\\.')
    deid <- (results %>% filter(padj<pval) %>% filter(abs(log2FoldChange)>l2fc))$transcript_id
    deset<- c(deset,deid)
    
    print(paste(control,condition,length(deid),length(deset)))
    
    # results<-results%>% left_join(annot, by = "gene_id")
    results %>% write.csv(paste0(res_dir_out,condition,'_vs_',control,'/DEseq_all_res.csv'), quote=F, row.names = F)
    
    #creating tables for GO
    go_lists <- results %>% filter(baseMean>0)
    rownames(go_lists) <- c()
    
    go_lists %>% select(transcript_id) %>%
      write.csv(paste0(res_dir_out,condition,'_vs_',control,'/background_list.txt'),row.names=F,col.names=F,quote=F)
    
    go_lists %>%
      filter(padj<pval) %>%
      filter(log2FoldChange>l2fc) %>%
      select(transcript_id) %>%
      write.csv(paste0(res_dir_out,condition,'_vs_',control,'/upregulated_list.txt'),row.names=F,col.names=F,quote=F)
    
    go_lists %>%
      filter(padj<pval) %>%
      filter(log2FoldChange<l2fc) %>%
      select(transcript_id) %>%
      write.csv(paste0(res_dir_out,condition,'_vs_',control,'/downregulated_list.txt'),row.names=F,col.names=F,quote=F)
    
    
    
    sigres <- results %>%
      filter(padj<pval) %>%
      filter(baseMean>0)%>%
      filter(abs(log2FoldChange)>l2fc)%>%
      arrange(desc(abs(log2FoldChange)),desc(baseMean))
    
    sigres %>% write.csv(file=paste0(res_dir_out,condition,'_vs_',control,'/diff_expression.csv'),  quote=F, row.names = F)
    
    #Numbers of differentially expressed genes
    n.differential = nrow(sigres)
    n.differential.up = sum(sigres$log2FoldChange > 0)
    n.differential.down = sum(sigres$log2FoldChange < 0)
    print(n.differential)
    print(n.differential.up)
    print(n.differential.down)
    
    
    grey = rgb(51, 51, 51,  50, maxColorValue = 255)
    red = rgb(255, 0, 0, 50, maxColorValue = 255)
    blue = rgb(0, 0, 255, 50, maxColorValue = 255)
    col.vec = rep(grey,nrow(res))
    col.vec[results$padj < pval & results$log2FoldChange < l2fc] = blue
    col.vec[results$padj < pval & results$log2FoldChange > l2fc] = red
    x_lim = max(abs(results$log2FoldChange), na.rm = T)
    pdf(paste0(res_dir_out,condition,'_vs_',control,'/Volcano_plot.pdf'))
    plot(results$log2FoldChange,-1*log10(results$padj),pch = 20,
         xlab = "Log2FoldChange", ylab = "-log10(padj)", col = col.vec, xlim = c(-1*x_lim, x_lim))
    legend("topright", paste(c("up = ","down = "),c(n.differential.up, n.differential.down), sep = ""), pch = 20, col = c(red,blue))
    dev.off()
    
  }
}



de_genes <- unique(deset)
write(de_genes, paste0(res_dir_out,'DE_at_any_timepoint.txt'))

#clustering DE genes ####
TPM <- read.csv('/project/Neurodifferentiation_System/Analysis_NGN3/Gene_expression_tables/GeneExpressionONT.csv', sep = ' ')
TPM <- TPM[-which(grepl("PAR_Y", TPM$gene, fixed = T)), ]

TPM<-TPM %>% separate(gene,into='gene_id',sep='\\.')
rownames(TPM) <- TPM$gene
TPM <- TPM[, sort(colnames(TPM))]
deTPM <- TPM[de_genes, ]
for (c in unique(df$condition)){
  deTPM[,paste0("Mean_",c)]<-deTPM%>% select(starts_with(c))%>%rowMeans()
}
meanTPM<-deTPM%>% select(starts_with("Mean"))
standTPM<-as.data.frame(t(scale(t(meanTPM))))
kclusters<-kmeans(standTPM,3)
standTPM$max<- sapply(as.data.frame(t(standTPM)),max)
standTPM$maxi<-max.col(standTPM)
standTPM$class<-kclusters$cluster
standTPM<- standTPM %>% arrange(desc(class))
kclasses<-standTPM %>%select(class)
standTPM<-standTPM%>%select(-c(class,max,maxi))
#standTPM<-standTPM%>%select(-c(class))

pheatmap(standTPM,show_rownames=F,cluster_cols=F,cluster_rows = F, annotation_colors=list("class"=qualitative_hcl(4,palette="Dark3")),
         annotation_row=kclasses, annotation_names_row = F, 
         labels_col = c('day0', 'day3', 'day5'), angle_col = '0',
         color =rev(scico(100,palette='roma'))[10:90], annotation_legend = F, fontsize_col = 20,
         filename = paste0(plot_dir_out, 'DE_heatmap.pdf'), width = 8, height = 5
         )




###How many novel isoforms are differentially expressed? ####
# fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter.1/nanopore_RulesFilter_result_classification.txt" 
fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_RulesFilter_result_classification.txt" 
df <- read.delim(fn)

filter_df <- df %>% filter(filter_result == "Isoform")
filter_df <- filter_df %>% filter(structural_category %in% c("full-splice_match", "incomplete-splice_match", "novel_in_catalog" ,"novel_not_in_catalog" ,"fusion"))

dte_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/DTE_analysis.1/DE_at_any_timepoint.txt"
# dte_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONTseq_quantification/Results/DTE_analysis/DE_at_any_timepoint.txt"
dte <- read.delim(dte_fn, header = F)

dte <- dte$V1

filter_df <- filter_df %>% mutate(DTE = ifelse(isoform %in% dte, "Yes", "No"))

sum(filter_df$DTE == "Yes")

sum(filter_df$DTE == "Yes" & (filter_df$structural_category == "novel_in_catalog" | filter_df$structural_category == "novel_not_in_catalog" | filter_df$structural_category == "fusion"))

sum(filter_df$structural_category == "novel_in_catalog" | filter_df$structural_category == "novel_not_in_catalog" | filter_df$structural_category == "fusion")
