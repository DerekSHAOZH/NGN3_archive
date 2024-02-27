library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(venn)
library(tximport)
library(patchwork)
library(DESeq2)
library(colorspace)
library(RColorBrewer)
library(pheatmap)
library(rtracklayer)
library(scico)
library(cluster)
library(factoextra)
plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Plots/'
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)


###QC: ONT data only ####
##load ONT gene count data ####
samples <- c('day0_1', 'day0_2', 'day0_3', 'day3_1', 'day3_2', 'day3_3', 'day5_1', 'day5_2', 'day5_3')
data_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts.txt'
data <- read.csv(file=data_fn, header=TRUE, sep=",")
data <- data[, sort(colnames(data))]
countData <- data[colnames(data)[2:(length(samples) + 1)]]

countData <- as.matrix(countData)

name <- colnames(data)[2:(length(samples) + 1)]
condition <- sapply(strsplit(name,"_"), `[`, 1)
colData <- cbind(name,condition)

##DESeq2 median-of-ratio norm and vst ####
dds <- DESeqDataSetFromMatrix(countData,colData,design=~condition)
dds <- estimateSizeFactors(dds)
Pcounts <- counts(dds, normalized=TRUE)

vsd <- vst(dds)
##1.heatmap of correlation ####
sampleDists <- dist(t(assay(vsd)))
sampleDists <- cor(assay(vsd))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- name
colnames(sampleDistMatrix) <- name
corrMatrix<-sampleDistMatrix


pheatmap(corrMatrix,col=colors, cluster_rows = F, cluster_cols = F,
         filename = paste0(plot_dir_out,'long_pearson_corr.pdf'))

##2.PCA ####
pdf(paste0(plot_dir_out,'long_PCA_count.pdf'), 6, 4)
# g <- plotPCA(vsd) +
#   theme_bw()
# print(g)
plotPCA(vsd) +
    theme_bw()
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  theme_bw() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  geom_point(size=4) +
  # labs(title = "Principle Component Analysis") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
        axis.text = element_text(face = "bold", size=12), axis.title = element_text(face = "bold", size=13),
        legend.title = element_text(size = 13, face = "bold"), legend.text = element_text(size = 12))


# ggsave(filename = paste0(plot_dir_out,'long_PCA.pdf'), plot = g,
#        device = 'pdf', width = 6, height = 6)
dev.off()


##3.switch per gene ####
#load ONT background gene list
# ref_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts.txt'
# ref <- read.csv(file=ref_fn, header=TRUE, sep=",")
# ref <- ref[, c((ncol(ref)-4):ncol(ref))]
# annot_ref <- ref[-which(ref$gene_id == ""), ]

long_ref_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/GffCompare/nanopore.combined.filt.gtf'
gtf <- rtracklayer::import(long_ref_fn)
gtf_df=as.data.frame(gtf)
gtf_df <- gtf_df %>% filter(type == 'transcript')


#annotated: with cmp_ref & class code == "="
annot_ref <- gtf_df[!is.na(gtf_df$cmp_ref) & gtf_df$class_code == '=',]
annot_ref <- annot_ref %>% select(transcript_id, gene_name, cmp_ref)

#total number of genes used: 10950
length(unique(annot_ref$gene_name))

#total number of isoforms used:18345
length(unique(annot_ref$cmp_ref))

days<-c(0,3,5)

ctrl <- sapply(days[1:2], function(x){y <- paste0("day", x)})
cond <- sapply(days[2:3], function(x){y <- paste0("day", x)})

#load ONT switch summary
data_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/switch_summary.csv"
data <- read.csv(data_fn)

#filter only subsequent switches
data <- data[data$condition_1 == ctrl[1] & data$condition_2 == cond[1] | 
               data$condition_1 == ctrl[2] & data$condition_2 == cond[2], ]

#filter significant isoform switch
data <- data[abs(data$dIF) > 0.1 & data$isoform_switch_q_value < 0.05, ]
data <- data %>% left_join(annot_ref, c('isoform_id' = 'transcript_id'))

#keep only annotated genes
data <- data[-which(is.na(data$gene_name.y)), ]

#number of genes with significant isoform switch: 295
length(unique(data$gene_name.y))

#number of isoform with significant switch events: 333
length(unique(data$isoform_id))

#number of significant isoform switch events: 344
length(data$isoform_id)

gene <- unique(data$gene_name.y)

#class of number of switch per gene
count <- c()
for (i in 1:length(gene)) {
  tmp <- sum(data$gene_name.y == gene[i])
  if(tmp == 4) {print(gene[i])}
  if (!(tmp %in% count)) {
    count <- c(count, tmp)
  }
}
count <- c(0, count)
count <- sort(count)

#number of genes in each class
df <- data.frame(num_switch = count, num_gene = c(rep(0, length(count))))
df[df$num_switch == 0, ]$num_gene = length(unique(annot_ref$gene_name)) - length(gene)
for (i in 1:length(gene)) {
  tmp <- sum(data$gene_name.y == gene[i])
  df[df$num_switch == tmp, ]$num_gene = df[df$num_switch == tmp, ]$num_gene + 1
}


#make barplot with gapped y axis
#whole figure
ggplot(data=df,aes(as.factor(num_switch),num_gene,fill=as.factor(num_switch)))+
  geom_bar(stat = "identity")+
  theme_classic()


#first fig
p1=ggplot(data=df,aes(as.factor(num_switch),num_gene,fill=as.factor(num_switch)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=num_gene), vjust=0, size = 5) + 
  coord_cartesian(ylim = c(0,300))+ #setting the down part
  theme_classic() +
  theme(legend.position="none",axis.title=element_text(size=16)) +
  xlab("Number of differential isoform expression events per gene") +
  ylab("Number of annotated gene")
  # ylab("Number of DE gene")


#second fig
p2=ggplot(data=df,aes(as.factor(num_switch),num_gene,fill=as.factor(num_switch)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=num_gene), vjust=0, size = 5) + 
  coord_cartesian(ylim = c(10000,11000))+ #setting the upper part
  guides(x = "none")+ # remove x line
  # labs(title = "Graph with broken y axis")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "none", # remove legend
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

#patchwork
#bind fig1+fig2
pdf(paste0(plot_dir_out, 'long_isoswitch_per_gene.pdf'), 7, 5)
# pdf(paste0(plot_dir_out, 'isoswitch_per_DE_gene.pdf'), 8, 6)
p2/p1
dev.off()

##4.annotated and novel transcript ####
#gtf annotation file but lack ref gene_id
long_ref_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/GffCompare/nanopore.combined.filt.gtf'
gtf <- rtracklayer::import(long_ref_fn)
gtf_df=as.data.frame(gtf)
gtf_df <- gtf_df %>% filter(type == 'transcript')
# gtf_df <- gtf_df %>% select(transcript_id, gene_name, cmp_ref)

#annotated: with cmp_ref & class code == "="
#Distribution of two classes of transcripts ####
num_annot <- nrow(gtf_df[!is.na(gtf_df$cmp_ref) & gtf_df$class_code == '=',])
df <- data.frame(group = c("Annotated transcript", "Novel transcript"), 
                 value = c(num_annot, nrow(gtf_df)-num_annot))
df <- df %>% 
  mutate(perc = `value` / sum(`value`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

pdf(paste0(plot_dir_out, 'long_transcript_type.pdf'), 5, 3)
ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = "black") +
  geom_text(aes(label = value), color = c("white", 1),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE,size = 8
            ) +
  guides(fill = guide_legend(title = "ONT-derived transcript", size = 10)) +
  scale_fill_viridis_d() +
  coord_polar(theta = "y") + 
  theme_void() +
  theme(legend.title=element_text(size=15), 
        legend.text=element_text(size=15))
dev.off()

#Distribution of five classes of transcripts from SQANTI3 ####
class_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_RulesFilter_result_classification.txt'
class_df <- read.csv(class_fn, sep = '\t')


type <- c('full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'fusion')

df <- data.frame(group = c("Full Splice Match (FSM)", "Incomplete Splice Match (ISM)", 
                           "Novel in Catalog (NIC)", "Novel Not in Catalog (NNC)", 'Fusion'), 
                 value = c(nrow(class_df[class_df$structural_category == type[1] & class_df$filter_result == 'Isoform',]),
                           nrow(class_df[class_df$structural_category == type[2] & class_df$filter_result == 'Isoform',]), 
                           nrow(class_df[class_df$structural_category == type[3] & class_df$filter_result == 'Isoform',]), 
                           nrow(class_df[class_df$structural_category == type[4] & class_df$filter_result == 'Isoform',]), 
                           nrow(class_df[class_df$structural_category == type[5] & class_df$filter_result == 'Isoform',])))
df$group <- factor(df$group, levels = c("Full Splice Match (FSM)", "Incomplete Splice Match (ISM)", 
                                        "Novel in Catalog (NIC)", "Novel Not in Catalog (NNC)", 'Fusion'))
df <- df %>% 
  mutate(perc = `value` / sum(`value`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc, accuracy = 0.1))

df$display_value <- sapply(df$value, function(x){
  if (x > 1000){
    if(x%%1000 < 100){
      y <- paste0(floor(x/1000), ',0', x%%1000)
    }else{
      y <- paste0(floor(x/1000), ',', x%%1000)
    }
  }else{
    y <- x
  }
  
  y
})

df$display_labels <- c('0.4%', '4.5%', '9.8%', '29.5%', '55.8%')

cols = c("Full Splice Match (FSM)" = "#6baed6", "Incomplete Splice Match (ISM)" = "#0EA293", 
         "Novel in Catalog (NIC)" = "#fc8d59", "Novel Not in Catalog (NNC)" = "#ee6a50",
         'Fusion' = '#FF2171')  
  
pdf(paste0(plot_dir_out, 'long_transcript_type_perc_sqanti3_legend.pdf'), 9, 6)
ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(color = "white") +
  geom_text(aes(label = paste0(display_value, '\n(', display_labels, ')')), color = c("white"),
            position = position_stack(vjust = 0.5),
            show.legend = FALSE,size = 7
  ) +
  # geom_text(aes(1.6, label = group), color = c("black"),
  #           position = position_stack(vjust = 0.5), hjust = "outward", 
  #           show.legend = FALSE,size = 7
  # ) +
  guides(fill = guide_legend(title = "Transcript isoform types", size = 10)) +
  scale_fill_manual(values = cols) +
  coord_polar(theta = "y", clip = "off") + 
  theme_void() +
  # theme(legend.position = "none")

  theme(legend.title=element_text(size=20),
        legend.text=element_text(size=20))
dev.off()

#robust novel transcripts ####

#join gtf with day0_1 coverage
cov <- read.table("/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/transcriptome_coverage_summary.txt", sep = "\t", header = T)
cov <- cov %>% rowwise() %>% mutate(mean_read = mean(c(day0_1,day0_2,day0_3,day3_1,day3_2,day3_3,day5_1,day5_2,day5_3)), 
                                    median_read = median(c(day0_1,day0_2,day0_3,day3_1,day3_2,day3_3,day5_1,day5_2,day5_3)))
gtf_read <- gtf_df %>% left_join(cov, by = c('oId' = 'rname') )
novel_gtf_read <- gtf_read[gtf_read$class_code != '=',]

###QC: comparison between short-read and ONT gene quantification ####
#load ONT gene count data ####
long_samples <- c('day0_1', 'day0_2', 'day0_3', 'day3_1', 'day3_2', 'day3_3', 'day5_1', 'day5_2', 'day5_3')
long_data_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts.txt'
long_data <- read.csv(file=data_fn, header=TRUE, sep=",")
long_data <- long_data[colnames(long_data)[c(2:(length(samples) + 1), which(colnames(long_data) == 'gene_id'))]]
long_data <- long_data[, sort(colnames(long_data))]



long_count <- aggregate(. ~ gene_id, data = long_data, FUN = sum)
long_count <- long_count[-1, ]
rownames(long_count) <- long_count$gene_id
long_count <- long_count[, -1]
colnames(long_count) <- unlist(lapply(colnames(long_count)[1:length(colnames(long_count))], function(x) {y <- paste0('long_', x)}))


#load short-read gene count data:  ####

dire<-'/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/'

samples <- c('OJ63','OJ64','OJ65','OJ66','OJ67','OJ68','OJ69','OJ70','OJ71','OJ72','OJ73','OJ74',
             'OJ75','OJ76','OJ77','OJ79','OJ80','OJ81','OJ83','OJ85','OJ86','OJ87','OJ88','OJ89')
cond <- c('day0','day0','day0','day0','day0','day1','day1','day1','day1','day1','day2','day2',
          'day2','day2','day2','day3','day3','day3','day3','day5','day5','day5','day5','day5')
rep <- c('rep1','rep2','rep3','rep4','rep5','rep1','rep2','rep3','rep4','rep5','rep1','rep2',
         'rep3','rep4','rep5','rep1','rep2','rep3','rep5','rep1','rep2','rep4','rep5','rep6')
fns <- paste0(dire,samples,'.genes.results')


days <- c(0, 3, 5)
index <- c()
for (i in 1:3) {
  index <- c(index, which(grepl(days[i], cond, fixed = T)))
}

samples_short <- samples[index]
cond_short <- cond[index]
rep_short <- rep[index]
short_fns <- paste0(dire,samples_short,'.genes.results')

short_data<-tximport(short_fns, type="rsem", txIn = FALSE, txOut = FALSE)
short_count <- short_data$counts
mode(short_count) <- "integer"
colnames(short_count) <- paste0('short_', cond_short, '_', rep_short)

#uniform gene id #######
short_count <- short_count[-which(grepl("PAR_Y", rownames(short_count), fixed = T)), ]
long_count <- long_count[-which(grepl("PAR_Y", rownames(long_count), fixed = T)), ]


for (i in 1:nrow(short_count)) {
  rownames(short_count)[i] <- unlist(strsplit(rownames(short_count)[i], '\\.'))[1]
}
for (i in 1:nrow(long_count)) {
  rownames(long_count)[i] <- unlist(strsplit(rownames(long_count)[i], '\\.'))[1]
}


#number of common genes: 14734
common <- intersect(rownames(short_count), rownames(long_count))
length(common)

short_count <- as.data.frame(short_count)
long_count <- as.data.frame(long_count)

short_count <- short_count[common, ]
long_count <- long_count[common, ]


count <- cbind(short_count, long_count)

count <- as.matrix(count)

#DESeq2 median-of-ratio normalization and vst ####
sample <- colnames(count)
condition <- c(paste0('short_', cond_short), paste0('long_', sapply(strsplit(long_samples,"_"), `[`, 1)))

colData <- cbind(sample,condition)
dds <- DESeqDataSetFromMatrix(count,colData,design=~condition)
dds <- estimateSizeFactors(dds)
Pcounts <- counts(dds, normalized=TRUE)

vsd <- vst(dds)

#1.heatmap of correlation ####
# sampleDists <- dist(t(assay(vsd)))
sampleDists <- cor(assay(vsd), method = "spearman")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- sample
colnames(sampleDistMatrix) <- sample
corrMatrix<-sampleDistMatrix
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
#pdf(paste0(plot_dir_out,'long_short_correlation.pdf'))
pheatmap(corrMatrix,col=colors, cluster_rows = F, cluster_cols = F,
         filename = paste0(plot_dir_out,'long_short_correlation.pdf')
         )
dev.off()

#2.PCA ####
pdf(paste0(plot_dir_out,'long_short_PCA.pdf'))
plotPCA(vsd) +
  theme_bw()
dev.off()






###QC: comparison between short-read and ONT DE genes ####
days<-c(0,3,5)

ctrl <- sapply(days[1:2], function(x){y <- paste0("day", x)})
cond <- sapply(days[2:3], function(x){y <- paste0("day", x)})

#stacked bar plot of three classes {short-only, ONT-only, both-short-and-ONT} ####

df <- data.frame(comparison = c(rep("day0 vs day3", 3), rep("day3 vs day5", 3)), 
                 group = rep(c("Both", "only in RNA-seq", "Only in ONT-seq"), 2),
                 value = 0)
for (i in 1:length(ctrl)) {
  long_DE_fn <- paste0('/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/', cond[i], '_vs_', ctrl[i], '/diff_expression.csv')
  long_DE_df <- read.csv(long_DE_fn, row.names = 1)
  
  short_DE_fn <- paste0('/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/', cond[i], '_background_', ctrl[i], '/diff_expression.csv')
  short_DE_df <- read.csv(short_DE_fn, row.names = 1)
  
  both <- intersect(short_DE_df$ensgene, long_DE_df$gene_id)
  short_only <- setdiff(short_DE_df$ensgene, both)
  long_only <- setdiff(long_DE_df$gene_id, both)
  
  df[df$comparison == paste0(ctrl[i], ' vs ', cond[i]), ]$value <- c(length(both), length(short_only), length(long_only))
  
}

pdf(paste0(plot_dir_out, 'bar_long_short_DE.pdf'), 8, 6)
ggplot(df, aes(x = comparison, y= value, fill = group)) +
  geom_bar(position = 'stack', stat = 'identity') +
  theme_bw()+
  xlab("Subsequent state transition") +
  ylab("Frequency of DE genes") +
  guides(fill=guide_legend(title="Occurrence"))
dev.off()


#venn diagram of DE genes between subsequent states ####
genes <- c(list(short_DE_df$ensgene), list(long_DE_df$gene_id))
names(genes) <- c('RNA-seq DE genes', 'ONT-seq DE genes')
pdf(paste0(plot_dir_out,'venn_long_short_DE.pdf'),6,4)
venn(genes, ilab=TRUE, zcolor = "style", box = F, ilcs = 1.3, sncs = 1)
dev.off()
###QC: comparison between short-read and ONT isoform switch events  ####
#load ONT isoformswitchAnalyzeR result ####
long_data_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/switch_summary.csv'
long_res <- read.csv(file=long_data_fn, header = T, sep = ',')

#gtf annotation file but lack ref gene_id
# long_ref_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/GffCompare/nanopore.combined.filt.gtf'
# gtf <- rtracklayer::import(long_ref_fn)
# gtf_df=as.data.frame(gtf)
# gtf_df <- gtf_df %>% filter(type == 'transcript')
# gtf_df <- gtf_df %>% select(transcript_id, gene_name, cmp_ref)

long_ref_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts.txt'
long_ref <- read.csv(file=long_ref_fn, header = T, sep = ',')
long_ref <- long_ref[!duplicated(long_ref$transcript_id), ]

diff <- setdiff(gtf_df$transcript_id, long_ref$transcript_id)

## LOOK INTO TCONS_00000001 and TCONS_00000002 in IGV!!! ####
## WHY 911 TRANSCRIPTS ARE MISSING IN GENE COUNT DATA!!! ####

# long_res <- long_res %>% left_join(gtf_df, by = c('isoform_id' = 'transcript_id'))
#load short-read isoformswitchAnalyzeR result
#common

#scatter plot of dIF ####
#venn diagram of genes with significant isoform switches ####
#stacked bar plot of three classes {short-only, ONT-only, both-short-and-ONT} ####
#venn diagram of significant isoform switches events ####
#stacked bar plot of three classes {short-only, ONT-only, both-short-and-ONT} ####

###QC: comparison of active genes/isoforms per day with short-read and ONT ####
##"Active" for certain day: TPM > 1 per replicate, then take median number  ####
##1. RNA-seq gene ####
#load short-read gene TPM data ####
gene_fn='/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_gene_TPM_shortreads_byDay.csv'
gene_df = read.csv(gene_fn)
#boxplot ####
df = data.frame(day = c(rep('0', 5), rep('1', 5), rep('2',5), rep('3',4), rep('5',5)),
                num = rep(0,24))
days = colnames(gene_df)[2:25]
for (i in 1:length(days)){
  # df$num[i] = sum(gene_df[days[i]] > 1)
  df$num[i] = sum(gene_df[days[i]] > 100)
}

# pdf(paste0(plot_dir_out, 'short_active_gene_1.pdf'), 5, 5)
pdf(paste0(plot_dir_out, 'short_active_gene_100.pdf'), 5, 5)
ggplot(df, aes(x = day, y = num))+
  geom_boxplot() +
  stat_compare_means() + 
  xlab("Day") +
  # ylab("RNA-seq: Number of active genes (TPM > 1)") +
  ylab("RNA-seq: Number of active genes (TPM > 100)") +
  theme_bw()
dev.off()

##2. RNA-seq isoform ####
iso_fn='/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_isoform_TPM_shortreads_byDay.csv'
iso_df = read.csv(iso_fn)

#boxplot ####
df = data.frame(day = c(rep('0', 5), rep('1', 5), rep('2',5), rep('3',4), rep('5',5)),
                num = rep(0,24))
days = colnames(iso_df)[3:26]
for (i in 1:length(days)){
  # df$num[i] = sum(iso_df[days[i]] > 1)
  df$num[i] = sum(gene_df[days[i]] > 100)
}

# pdf(paste0(plot_dir_out, 'short_active_isoform_1.pdf'), 5, 5)
pdf(paste0(plot_dir_out, 'short_active_isoform_100.pdf'), 5, 5)
ggplot(df, aes(x = day, y = num))+
  geom_boxplot() +
  stat_compare_means() + 
  xlab("Day") +
  # ylab("RNA-seq: Number of active isoforms (TPM > 1)") +
  ylab("RNA-seq: Number of active isoforms (TPM > 100)") +
  theme_bw()
dev.off()

##3. ONT-seq isoform ####
iso_fn ='/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_deseq2norm.txt'
iso_df = read.csv(iso_fn)
iso_df = iso_df[, order(names(iso_df))]
#boxplot ####
df = data.frame(day = c(rep('0', 3), rep('3',3), rep('5',3)),
                num = rep(0,9))
days = colnames(iso_df)[2:10]
for (i in 1:length(days)){
  # df$num[i] = sum(iso_df[days[i]] > 1)
  df$num[i] = sum(iso_df[days[i]] > 100)
}

# pdf(paste0(plot_dir_out, 'long_active_isoform_1.pdf'), 5, 5)
pdf(paste0(plot_dir_out, 'long_active_isoform_100.pdf'), 5, 5)
ggplot(df, aes(x = day, y = num))+
  geom_boxplot() +
  stat_compare_means() + 
  xlab("Day") +
  # ylab("ONT-seq: Number of active isoforms (TPM > 1)") +
  ylab("ONT-seq: Number of active isoforms (TPM > 100)") +
  theme_bw()
dev.off()

##3. ONT-seq gene ####
gene_fn ='/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_byGene_deseq2norm.csv'
gene_df = read.csv(gene_fn)

#boxplot ####
df = data.frame(day = c(rep('0', 3), rep('3',3), rep('5',3)),
                num = rep(0,9))
days = colnames(gene_df)[3:11]
for (i in 1:length(days)){
  # df$num[i] = sum(gene_df[days[i]] > 1)
  df$num[i] = sum(gene_df[days[i]] > 100)
}

# pdf(paste0(plot_dir_out, 'long_active_gene_1.pdf'), 5, 5)
pdf(paste0(plot_dir_out, 'long_active_gene_100.pdf'), 5, 5)
ggplot(df, aes(x = day, y = num))+
  geom_boxplot() +
  stat_compare_means() + 
  xlab("Day") +
  # ylab("ONT-seq: Number of active genes (TPM > 1)") +
  ylab("ONT-seq: Number of active genes (TPM > 100)") +
  theme_bw()
dev.off()

###QC: example of AS consequences ####
#load switch-wise consequence summary ####
conseq_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/consequence_detailed_sqanti.csv'
conseq_df <- read.csv(conseq_fn)

conseq_df <- conseq_df %>% filter(isoformsDifferent == TRUE)
conseq_df <- conseq_df %>% filter(featureCompared == 'ORF_seq_similarity' |
                                    featureCompared == 'domains_identified' |
                                    featureCompared == 'NMD_status' | 
                                    featureCompared == 'domain_length' |
                                    featureCompared == 'signal_peptide_identified')
conseq_df <- conseq_df %>% filter(Up_structural_category == 'full-splice_match' & Down_structural_category == 'full-splice_match')
conseq_df <- conseq_df %>% filter(Up_associated_gene == Down_associated_gene)


# test <- conseq_df %>% mutate(switch  = paste0(Comparison, '_', isoformUpregulated, '_', isoformDownregulated))
# test <- test[!duplicated(test$switch), ]
# sum(test$Comparison == 'day0 vs day3')
# sum(test$Comparison == 'day3 vs day5')


# conseq_df <- conseq_df[, 2:8]
# sq_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_classification_TPM.txt'
# sq_df <- read.csv(sq_fn, sep = '\t')
# df <- sq_df[, c(1,6,7,8)]
# 
# test <- conseq_df %>% left_join(df, by = c("isoformUpregulated" = "isoform"))
# test <- test %>% left_join(df, by = c("isoformDownregulated" = "isoform"))
# colnames(test)[which(colnames(test) == "structural_category.x")] <- "Up_structural_category"
# colnames(test)[which(colnames(test) == "associated_gene.x")] <- "Up_associated_gene"
# colnames(test)[which(colnames(test) == "associated_transcript.x")] <- "Up_associated_transcript"
# colnames(test)[which(colnames(test) == "structural_category.y")] <- "Down_structural_category"
# colnames(test)[which(colnames(test) == "associated_gene.y")] <- "Down_associated_gene"
# colnames(test)[which(colnames(test) == "associated_transcript.y")] <- "Down_associated_transcript"
# 
# write.csv(test, '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/consequence_detailed_sqanti.csv', quote = F)

#load ONT TPM data ####
count_fn ='/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_deseq2norm.txt'
count_df = read.csv(count_fn)
count_df = count_df[, order(names(count_df))]

###QC: example of SQANTI transcript types ####
sq_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_classification_TPM.txt'
sq_df <- read.csv(sq_fn, sep = '\t')
df <- sq_df[, c(1,6,7,8)]

test <- conseq_df %>% left_join(df, by = c("isoformUpregulated" = "isoform"))
test <- test %>% left_join(df, by = c("isoformDownregulated" = "isoform"))
colnames(test)[which(colnames(test) == "structural_category.x")] <- "UP_structure_category"
colnames(test)[which(colnames(test) == "associated_gene.x")] <- "Up_associated_gene"
colnames(test)[which(colnames(test) == "associated_transcript.x")] <- "Up_associated_transcript"
colnames(test)[which(colnames(test) == "structural_category.y")] <- "Down_structure_category"
colnames(test)[which(colnames(test) == "associated_gene.y")] <- "Down_associated_gene"
colnames(test)[which(colnames(test) == "associated_transcript.y")] <- "Down_associated_transcript"

test <- sq_df %>% filter(structural_category == 'fusion')
test <- test[order(test$FL_TPM, decreasing = T), ]


###QC:read covering exon-exon junctions in ONT-seq and RNA-seq ####
##PFN2 ####
summary_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/perc_read_cover_junction.PFN2.csv'
summary <- read.csv(summary_fn)
summary <- summary[, -1]

pdf(paste0(plot_dir_out, 'perc_read_cover_junction.PFN2.pdf'), 6, 5)
ggplot(summary, aes(x=num, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(.9)) +
  scale_fill_discrete(name="Read length",
                     labels=c("Long read (ONT-seq)","Short read (RNA-seq)")
                    ) + 
  labs(x = 'Number of exon-exon junctions',
       y = 'Percentage of mappable reads',
       title = 'PFN2') + 
  theme_bw() +
  theme(text = element_text(size = 15)) 
dev.off()
##whole genome ####
summary_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/perc_read_cover_junction.WG.csv'
summary <- read.csv(summary_fn)
summary <- summary[, -1]

pdf(paste0(plot_dir_out, 'perc_read_cover_junction.WG.pdf'), 7, 5)
ggplot(summary, aes(x=num, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(.9)) +
  scale_fill_discrete(name="Read length",
                      labels=c("Long read (ONT-seq)","Short read (RNA-seq)")
  ) + 
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  labs(x = 'Number of exon-exon junctions',
       y = 'Percentage of mappable reads',
       title = 'Whole genome') + 
  theme_bw() +
  theme(text = element_text(size = 15)) 
dev.off()

##only genes ####
summary_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/perc_read_cover_junction.gene.csv'
summary <- read.csv(summary_fn)
summary <- summary[, -1]

#remove reads cover 0 junction
summary <- summary %>% filter(num != 0)

summary$mean <- summary$mean * 100
summary$sem <- summary$sem * 100

pdf(paste0(plot_dir_out, 'perc_read_cover_junction.gene.pdf'), 7, 4)
ggplot(summary, aes(x=num, y=mean, fill=group)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem), width=.2,
                position=position_dodge(.9)) +
  scale_fill_discrete(name="Read length",
                      labels=c("Long read (ONT-seq)","Short read (RNA-seq)")
  ) + 
  scale_x_continuous(breaks=seq(1, 10, 1)) +
  labs(x = 'Number of exon-exon junctions',
       y = 'Uniquely mapped reads (%)',
       title = 'Multi-exonic protein-coding genes') + 
  theme_bw() +
  theme(text = element_text(size = 15),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) 
dev.off()

###Heat,ap of transcript isoform expression ####
#load ONT TPM data ####
count_fn ='/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_TPM.txt'
count_df = read.csv(count_fn)
count_df = count_df[, order(names(count_df))]
count_df = count_df[!duplicated(count_df$transcript_id), ]
#only valid isoform after filtering ####
valid_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_inclusion-list.txt'
valid_df <- read.csv(valid_fn, header = F)
valid <- valid_df[, 1]
count_df <- count_df[count_df$transcript_id %in% valid, ]

rownames(count_df) <- count_df$transcript_id

day <- c('day0', 'day3', 'day5')
for (c in day){
  count_df[,paste0("Mean_",c)]<-count_df%>% select(starts_with(c))%>%rowMeans()
}


meanTPM<-count_df%>% select(starts_with("Mean"))

##z-scale the mean TPM for each transcript ####
standTPM<-as.data.frame(t(scale(t(meanTPM))))

##find optimal k for kmeans clustering ####
silhouette_score <- function(k){
  km <- kmeans(standTPM, centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(standTPM))
  mean(ss[, 3])
}
k <- 2:6
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)

fviz_nbclust(standTPM, kmeans, method = "wss")
fviz_nbclust(standTPM, kmeans, method = "silhouette")
fviz_nbclust(standTPM, kmeans, method = "gap_stat")


kclusters<-kmeans(standTPM,centers = 3, nstart=25)
standTPM$max<- sapply(as.data.frame(t(standTPM)),max)
standTPM$maxi<-max.col(standTPM)
standTPM$class<-kclusters$cluster
standTPM<- standTPM %>% arrange(desc(class))


class_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_RulesFilter_result_classification.txt'
class_df <- read.csv(class_fn, sep = '\t')

kclasses <- read.csv('/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/cluster_transcript.only_valid.csv',
                     row.names = 1)

# kclasses<-standTPM %>%select(class)
kclasses$class <- factor(kclasses$class, levels = c(1,2,3))
kclasses$transcript_id <- rownames(kclasses)
kclasses <- kclasses %>% left_join(class_df[, c('isoform', 'structural_category', 'filter_result')], 
                                   by = c('transcript_id' = 'isoform'))
type <- c('full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'fusion')
type_new <- c("Full Splice Match", "Incomplete Splice Match", 
              "Novel in Catalog", "Novel Not in Catalog", 'Fusion')
for (i in 1:length(type)) {
  kclasses[kclasses$structural_category == type[i], ]$structural_category <- type_new[i]
}
kclasses[!(kclasses$structural_category %in% type_new), ]$structural_category <- 'Other'

kclasses$filter_result <- factor(kclasses$filter_result, levels = c("Isoform", "Artifact"))
kclasses$structural_category <- factor(kclasses$structural_category, levels = c(type_new,'Other'))
# kclasses$class <- factor(kclasses$class, levels = rev(c(3,2,1)))
rownames(kclasses) <- kclasses$transcript_id
kclasses <- kclasses[, -which(colnames(kclasses) == 'transcript_id')]
standTPM<-standTPM%>%select(-c(class,max,maxi))
#standTPM<-standTPM%>%select(-c(class))


old_1 <- which(kclasses$class == 1)
old_2 <- which(kclasses$class == 2)
old_3 <- which(kclasses$class == 3)
kclasses$class[old_2] <- 1
kclasses$class[old_1] <- 2
kclasses$class <- factor(kclasses$class, levels = c(1,2,3))


tmp <- with(kclasses, kclasses[order(class, structural_category, decreasing = T), ]) 
tmp <- kclasses %>% arrange(factor(class, levels = c(3,2,1)), 
                            factor(filter_result, levels = c('Isoform', 'Artifact')),
                            factor(structural_category, levels = c(type_new,'Other')))
test <- standTPM[rownames(tmp ), ]


##heatmap ####

Var1        <- sequential_hcl(n=3, palette = 'Batlow')
names(Var1) <- factor(c(1,2,3))
Var2 <- c("Full Splice Match" = "#6baed6", "Incomplete Splice Match" = "#0EA293", 
                 "Novel in Catalog" = "#fc8d59", "Novel Not in Catalog" = "#ee6a50",
                 'Fusion' = '#FF2171', 'Other' = '#AEC3AE')
Var3 <- c('Artifact'= '#FF7D7D', 'Isoform'= '#59CE8F')
anno_colors <- list('class' = Var1,
                    'structural_category' = Var2,
                    'filter_result' = Var3)


pheatmap(test,show_rownames=F,cluster_cols=F,cluster_rows = F, annotation_colors=anno_colors,
         annotation_row=kclasses, annotation_names_row = F, 
         labels_col = c('day 0', 'day 3', 'day 5'), angle_col = '0',
         color =rev(scico(100,palette='roma'))[10:90],
         annotation_legend = T, fontsize_col = 12, fontsize_row = 5, border_color = NA, 
         # filename = paste0(plot_dir_out, 'transcript_heatmap_ONTseq.only_valid.pdf'), width = 8, height = 5
)

write.csv(kclasses, '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/cluster_transcript.only_valid.csv',
          quote = F)

