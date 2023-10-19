list.of.packages <- c("data.table","ggplot2","ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
MS.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_noSingleExon/de-novo/IsoformIntensities.txt"

args <- commandArgs(trailingOnly = TRUE) 

MS.path = args[1]
out = args[2]
out2 = args[3]
out3=args[4]

new_MS = as.data.table(read.table(MS.path, header = T, sep="\t"))

############################# median for replicate measurements
collapsed_MS = new_MS[,.(meanLFQ=median(meanLFQ),
                         sumLFQ=median(sumLFQ),
                         mean_norm_1=median(mean_norm_1),
                         sum_norm_1=median(sum_norm_1),
                         gene_name = unique(gene_name),
                        gene_id = unique(gene_id)),
                      by=c("transcript_id","Day")]

conditions=c("Day0","Day1","Day2","Day3","Day4","Day5") 
logFC_condition=c("Day1vsDay0","Day2vsDay1","Day3vsDay2","Day4vsDay3","Day5vsDay4")
# add NA for missing values
logFC_matrix=collapsed_MS[,.(Day=c(Day, conditions[!(conditions %in% Day)]), 
                             mean_norm_1=c(mean_norm_1,rep(NA,length(conditions)-length(mean_norm_1))), 
                             sum_norm_1=c(sum_norm_1,rep(NA,length(conditions)-length(sum_norm_1)))), 
                             by=c("transcript_id","gene_id", "gene_name")]

logFC_matrix=logFC_matrix[,.(mean_Day1vsDay0=log2(mean_norm_1[Day=="Day1"]/mean_norm_1[Day=="Day0"]),
                             mean_Day2vsDay1=log2(mean_norm_1[Day=="Day2"]/mean_norm_1[Day=="Day1"]),
                             mean_Day3vsDay2=log2(mean_norm_1[Day=="Day3"]/mean_norm_1[Day=="Day2"]),
                             mean_Day4vsDay3=log2(mean_norm_1[Day=="Day4"]/mean_norm_1[Day=="Day3"]),
                             mean_Day5vsDay4=log2(mean_norm_1[Day=="Day5"]/mean_norm_1[Day=="Day4"]),
                             sum_Day1vsDay0=log2(sum_norm_1[Day=="Day1"]/sum_norm_1[Day=="Day0"]),
                             sum_Day2vsDay1=log2(sum_norm_1[Day=="Day2"]/sum_norm_1[Day=="Day1"]),
                             sum_Day3vsDay2=log2(sum_norm_1[Day=="Day3"]/sum_norm_1[Day=="Day2"]),
                             sum_Day4vsDay3=log2(sum_norm_1[Day=="Day4"]/sum_norm_1[Day=="Day3"]),
                             sum_Day5vsDay4=log2(sum_norm_1[Day=="Day5"]/sum_norm_1[Day=="Day4"])), 
                          by=c("transcript_id","gene_id", "gene_name")]

logFC_table=logFC_matrix[,.(Day=logFC_condition,
                            mean_logFC=c(mean_Day1vsDay0, mean_Day2vsDay1, mean_Day3vsDay2, mean_Day4vsDay3, mean_Day5vsDay4),
                            sum_logFC=c(sum_Day1vsDay0, sum_Day2vsDay1, sum_Day3vsDay2, sum_Day4vsDay3, sum_Day5vsDay4)),
                         by=c("transcript_id","gene_id", "gene_name")]

write.table(logFC_matrix, col.names = T, row.names = F, sep = "\t", quote = F, file=out)

pdf(out3, width=12)
for (gene in unique(logFC_table$gene_id)){
  selection=logFC_table[logFC_table$gene_id==gene,]
  gene_name=unique(selection$gene_name)
  
  p2=ggplot(data=selection, aes(x=Day, y=mean_logFC, fill =transcript_id)) +
    geom_hline(yintercept=0) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)
    ) +
    geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5) +
    ylim(-5, 5) +
    labs(title=paste(gene_name," (",gene,")",sep=""), y="Trend Change") +
    theme_bw() +
    theme(legend.title=element_blank(), axis.title.x=element_blank()) 
  
  
  p=ggarrange(p2, 
              nrow = 1, ncol=1,
              align="hv",
              common.legend = TRUE, legend="top")
  
  print(p)
}
dev.off()



pdf(out2,width=12)
for (gene in unique(logFC_table$gene_id)){
  selection=logFC_table[logFC_table$gene_id==gene,]
  gene_name=unique(selection$gene_name)
  
  p1=ggplot(data=selection, aes(x=Day, y=sum_logFC, fill =transcript_id)) +
    geom_hline(yintercept=0) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)
    ) +
    geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5) +
    ylim(-5, 5) +
    labs(title=paste(gene_name," (",gene,")",sep="")) +
    theme_bw() +
    theme(legend.title=element_blank(), axis.title.x=element_blank()) 
    
  p2=ggplot(data=selection, aes(x=Day, y=mean_logFC, fill =transcript_id)) +
    geom_hline(yintercept=0) +
    stat_summary(
      fun.y = median,
      size=2,
      geom = 'line',
      aes(group = transcript_id, colour = transcript_id)
    ) +
    geom_dotplot(binaxis='y', stackdir='center',dotsize=0.5) +
    ylim(-5, 5) +
    labs(title=paste(gene_name," (",gene,")",sep="")) +
    theme_bw() +
    theme(legend.title=element_blank(), axis.title.x=element_blank()) 
  
  
  p=ggarrange(p1,p2, 
              nrow = 1, ncol=2,
              align="hv",
              common.legend = TRUE, legend="top")
  
  print(p)
}
dev.off()


