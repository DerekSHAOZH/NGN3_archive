list.of.packages <- c("data.table","ggplot2","ggpubr","pheatmap","ggcorrplot","tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("ggcorrplot"))
suppressPackageStartupMessages(library("tidyr"))


ratio.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/MSratio.atLeastTwoProteinIsoforms.ont.tab"
overlap.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/overlap.tab"
RNA.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/ONT-seq.txt"
out.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/continuous/trendChanges.ONT."

args <- commandArgs(trailingOnly = TRUE) 


ratio.path = args[1]
overlap.path = args[2]
RNA.path = args[3]
out=args[4]

RATIO=fread(ratio.path, header=T)
overlap=read.table(overlap.path)


RNA=fread(RNA.path, header=T)[,c("log2FoldChange","transcript_id","Comparison")]
colnames(RNA)=c("RNA_log2FC","transcript_id","Comparison")

TEMP=RATIO[,.(transcript_id_split=unlist(strsplit(transcript_id,split=";")),
              gene_id=gene_id,
              gene_name=unique(gene_name),
              groupSize=length(unlist(strsplit(transcript_id,split=";")))),by="transcript_id"]

TEMP=merge(TEMP,overlap, by.x="transcript_id_split", by.y="V1")
TEMP=TEMP[TEMP$V2>TEMP$groupSize,]
TEMP=merge(TEMP,TEMP[,.(Count=.N),by="gene_id"], by="gene_id")
TEMP=subset(TEMP,Count>=2)

genes_of_interest=unique(TEMP$gene_id)
transcripts_of_interest=unique(TEMP$transcript_id)

RATIO_FILTER=subset(RATIO, transcript_id %in% transcripts_of_interest)

####
RATIO_FILTER_TIDY=gather(RATIO_FILTER, key=Comparison, value=MS_log2FC, colnames(RATIO_FILTER)[grepl("Day",colnames(RATIO_FILTER))])
RATIO_FILTER_TIDY=merge(RATIO_FILTER_TIDY,RNA,by=c("transcript_id","Comparison"), all.x=T)


myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

RATIO_FILTER=myspread(RATIO_FILTER_TIDY,Comparison, c(MS_log2FC,RNA_log2FC))
colnames(RATIO_FILTER)=gsub("_log2FC","",colnames(RATIO_FILTER))
RATIO_FILTER=RATIO_FILTER[ , c(grep("MS",colnames(RATIO_FILTER)), grep("RNA",colnames(RATIO_FILTER))) ]

print("n=")
nrow(RATIO_FILTER)
CORR=cor(as.data.frame(RATIO_FILTER)[,grepl("Day",colnames(RATIO_FILTER))], method = "spearman")

CORR=CORR[grep("RNA",row.names(CORR)),grep("MS",row.names(CORR))]
pdf(out)
ggcorrplot(CORR, method = "circle",lab = TRUE, digits=1, legend.title="Spearman\nCorrelation")
dev.off()
