list.of.packages <- c("data.table","DEP","dplyr","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("DEP"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("stringr"))

data.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/MS_NeuroDifferentation_NGN3/proteinGroups.csv"
annotation.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS_SQANTI3_customProteome/consecutive/annotation.info.txt"
comparisonL="Day1_vs_Day0;Day2_vs_Day1;Day3_vs_Day2;Day4_vs_Day3;Day5_vs_Day4"
mode="atLeastTwoProteinIsoforms"
selection="unique"
args <- commandArgs(trailingOnly = TRUE) 

data.path = args[1]
annotation.path = args[2]
comparisonL = args[3] 
selection=args[4]
mode =args[5]
out = args[6]

comparison=unlist(strsplit(comparisonL, split=";"))

data=read.table(data.path,sep="\t", header = T)
# get selection
TEMP2=as.data.table(data)[,.SD,.SDcols=c("Peptide.counts..razor.unique.",paste("Peptide.counts..",selection,".",sep=""))]
data=data[TEMP2[,1]==TEMP2[,2],]

#remove potential contaminats
data <- subset(data, Reverse!="+" & Potential.contaminant!="+")

LFQ_columns <- grep("LFQ.", colnames(data)) # get LFQ column numbers

data=data[,c(1,LFQ_columns)]
annotation=read.table(annotation.path, header =T)


#add gene id
TEMP=as.data.table(data)
TEMP=TEMP[,.(Protein.ID=unlist(strsplit(Protein.IDs,split=";"))),by="Protein.IDs"]

TEMP=merge(TEMP, annotation, by.y="transcript_id", by.x="Protein.ID")

TEMP=TEMP[,.(Gene.IDs=paste(unique(gene_id),collapse=";"),
             Gene.names=paste(unique(gene_name),collapse=";")),by="Protein.IDs"]

data=unique(merge(data,TEMP,by="Protein.IDs"))

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.IDs", "Protein.IDs", delim = ";")

#TODO: adapt experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers

experimental_design <- data.frame(label=substring(colnames(data_unique)[LFQ_columns],15),
                                  condition=do.call(rbind,strsplit(colnames(data_unique)[LFQ_columns],split="\\."))[,3],
                                  replicate=rep(1:3,times=6))
data_se <- make_se(data_unique, LFQ_columns, experimental_design)


#plot_frequency(data_se)
data_filt <- filter_missval(data_se, thr = 1)
#plot_numbers(data_filt)
#plot_coverage(data_filt)
data_norm <- normalize_vsn(data_filt)
#plot_normalization(data_filt, data_norm)
#plot_missval(data_filt)
#plot_detect(data_filt)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
#plot_imputation(data_norm, data_imp_man)

data_diff_manual <- test_diff(data_imp_man, type = "manual", 
                              test = comparison)

dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = 1)
#plot_cond(dep)

data_results <- get_results(dep)

ratio_columns <- grep("_ratio", colnames(data_results)) # get LFQ column numbers
data_results=merge(data_unique[,c("Protein.IDs","Gene.IDs","name","Gene.names")],data_results[,c(1,ratio_columns)],by="name")[,-1]

genes_of_interest=names(table(data_results$Gene.IDs)[table(data_results$Gene.IDs)>=2])

if(mode=="atLeastTwoProteinIsoforms"){
  data_results=subset(data_results,Gene.IDs %in% genes_of_interest)
}

write.table(data_results, row.names = F, col.names = T, sep = "\t", quote = F, file=out)

