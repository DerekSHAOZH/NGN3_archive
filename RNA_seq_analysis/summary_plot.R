library(dplyr)
library(tidyr)
library(ggplot2)
library(venn)
library(tximport)
library(patchwork)
library(rtracklayer)
library(ggpolypath)
library(purrr)
library(scico)
library(colorspace)
library(RColorBrewer)
library(pheatmap)
library(conflicted)
library(cluster)
library(factoextra)
library(DESeq2)
library(tidyverse)
library(ggprism)
library(rstatix)
library(ggpubr)
library(scales)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# all_dire <- "/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitchAnalyser/GenesWithSwitches/"
all_dire <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/"
# filt_dire <- "/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitchAnalyser/FilteredMain/"
# data_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitchAnalyser/Results/switchSummary.csv"
data_fn = "/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/switch_summary.csv"
DE_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/DE_at_any_timepoint.txt'
event_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/AS_type_summary.csv"

#ref_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/OJ63.isoforms.results"
#ref_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_isoform_TPM_shortreads_byDay.csv"
ref_fn <- '/project/owlmayer/Applications/06_IsoformPipeline/data/genomeGencode/gencode.v28.annotation.gtf'

plot_dir_out <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Plots/"
res_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/'


days<-c(0,1,2,3,5)
###0.Data preparation ####
ctrl <- sapply(days[1:4], function(x){y <- paste0("day", x)})
cond <- sapply(days[2:5], function(x){y <- paste0("day", x)})

#load short-read background gene list ####
# ref <- read.csv(ref_fn)
# ref <- ref[, c(1, 2, ncol(ref)-1, ncol(ref))]

#gtf annotation file 
gtf <- rtracklayer::import(ref_fn)
gtf_df=as.data.frame(gtf)
gtf_df <- gtf_df %>% filter(type == 'transcript')
gtf_df <- gtf_df[-which(grepl('PAR_Y', gtf_df$gene_id, fixed = T)), ]
ref <- gtf_df %>% dplyr::select(gene_id, gene_name, gene_type)
ref <- ref[!duplicated(ref$gene_name), ]

gene_list <- unique(ref$gene_id)

background_gene <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/isoformswitchanalyzer_prefilter_gene.csv", row.names = 1)
background_gene <- background_gene$x

background_ref <- ref%>% filter(gene_name %in% background_gene)
#total number of background genes used in isoformswitchanalyzer: 14655
length(unique(background_ref$gene_id))
# write.table(background_ref$gene_id, "/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/isoformswitchanalyzer_prefilter_geneid.txt",
#             quote = F, row.names = F, col.names = F )

#total number of genes used: 56832

#load short-read DE gene list ####
DE_gene <- read.csv(DE_fn, header = F)
colnames(DE_gene) <- 'gene_id'
DE_gene <- DE_gene[-which(grepl('PAR_Y', DE_gene$gene_id, fixed = T)), ]
DE_gene <- as.data.frame(DE_gene)
colnames(DE_gene) <- 'gene_id'
# DE_gene <- DE_gene %>% separate(gene_id, into='gene_id', sep = '\\.')

#filter only DE genes as background
DE_ref <- background_ref[background_ref$gene_id %in% DE_gene$gene_id, ]
NotDE_ref <- background_ref[!(background_ref$gene_id %in% DE_gene$gene_id), ]

#total number of DE genes: 19892
length(unique(DE_gene$gene_id))

#total number of DE genes used in isoformswitchanalyzer: 10979
length(unique(DE_ref$gene_id))

#load short-read switch summary ####
# data <- read.csv(data_fn, sep = "\t")
data <- read.csv(data_fn)

#day0 vs day3 for comparison with ONT-seq
# tmp <- data[data$condition_1 == 'day0' & data$condition_2 == 'day3', ]
# tmp <- tmp[abs(tmp$dIF) > 0.1 & tmp$isoform_switch_q_value < 0.05, ]
# 
# tmp <- tmp %>% left_join(ref, by = "gene_name")
# list <- tmp$gene_id.y
# list <- as.data.frame(list)
# write.csv(list, file = paste0(res_dir_out, 'GenesWithIsoSwitch_day0_vs_day3.csv'), row.names = FALSE, col.names = NA, quote = FALSE)

#filter only subsequent switches
# data <- data[data$condition_1 == ctrl[1] & data$condition_2 == cond[1] | 
#                data$condition_1 == ctrl[2] & data$condition_2 == cond[2] |
#                data$condition_1 == ctrl[3] & data$condition_2 == cond[3] |
#                data$condition_1 == ctrl[4] & data$condition_2 == cond[4], ]

#filter significant isoform switches
data <- data[abs(data$dIF) > 0.1 & data$isoform_switch_q_value < 0.05, ]

all_dat <- data %>% left_join(ref, by = "gene_name")



##genes with transcript usage changes in background genes ####
background_data <- data %>% left_join(background_ref, by = "gene_name")


# background_data <- background_data[!is.na(background_data$gene_id.y), ]
#number of genes with significant isoform switch from any timepoint: 9291
length(unique(background_data$gene_id.y))

#number of isoform with significant switch events: 23930
length(unique(background_data$isoform_id))

#number of significant isoform switch events: 69344
length(background_data$isoform_id)

background_gene_change <- unique(background_data$gene_id.y)




##genes with transcript usage changes in DE genes ####
DE_data <- data %>% left_join(DE_ref, by = "gene_name")


DE_data <- DE_data[!is.na(DE_data$gene_id.y), ]
#number of DE genes with significant isoform switch from any timepoint: 6966
length(unique(DE_data$gene_id.y))

#number of isoform with significant switch events: 17959
length(unique(DE_data$isoform_id))

#number of significant isoform switch events: 52056
length(DE_data$isoform_id)

DE_gene_change <- unique(DE_data$gene_id.y)


##genes with transcript usage changes in NotDE genes ####
NotDE_data <- data %>% left_join(NotDE_ref, by = "gene_name")


NotDE_data <- NotDE_data[!is.na(NotDE_data$gene_id.y), ]
#number of NotDE genes with significant isoform switch from any timepoint: 2325
length(unique(NotDE_data$gene_id.y))

#number of isoform with significant switch events: 5971
length(unique(NotDE_data$isoform_id))

#number of significant isoform switch events: 17288
length(NotDE_data$isoform_id)

NotDE_gene_change <- unique(NotDE_data$gene_id.y)

#output genes with significant subsequent switches ####
# for (i in 1:length(ctrl)) {
#   list <- data[data$condition_1 == ctrl[i] & data$condition_2 == cond[i], ]$gene_id.y
#   list <- as.data.frame(list)
#   write.csv(list, file = paste0(res_dir_out, 'GenesWithIsoSwitch_', ctrl[i], '_vs_', cond[i], '.csv'), row.names = FALSE, col.names = NA, quote = FALSE)
#   
# }



###1.Overlapping between DE genes and genes with switch(s) ####
genes <- c(list(unique(DE_gene$gene_id)), list(gene))
names(genes) <- c('DE genes', 'Genes with isoSwitch')
pdf(paste0(plot_dir_out,'venn_DE_switch.pdf'),6,4)
venn(genes, ilab=TRUE, zcolor = "style", box = F, ilcs = 1.3, sncs = 1)
dev.off()


dom_switch_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/dominant_switch_gene.csv'
dom_switch <- read.csv(dom_switch_fn, header = F)
DE_list<-as.vector(as.character(map(strsplit(DE_gene$gene_id,'\\.'),1)))
iso_switch_list <-as.vector(as.character(map(strsplit(gene,'\\.'),1)))
genes <- c(list(DE_list), list(dom_switch$V1), list(iso_switch_list))
names(genes) <- c('DE genes', 'Genes with domIsoChange', 'Genes with isoSwitch')
pdf(paste0(plot_dir_out,'venn_DE_switch_domIsoChange.pdf'),6,4)
venn(genes, ilab=TRUE, zcolor = "style", box = F, ilcs = 1, sncs = 0.85)
dev.off()

tmp <- setdiff(intersect(dom_switch$V1, iso_switch_list), DE_list)
# write.csv(tmp, file = paste0(res_dir_out, 'mainIsoNotDE.csv'), row.names = FALSE, col.names = NA, quote = FALSE)


# gene <- intersect(unique(ref$gene_name), gene)

###2.Switch per gene ####
#find class of number of switch per gene ####
count <- c()

for (i in 1:length(gene)) {
  tmp <- sum(data$gene_id.y == gene[i])
  # if(tmp == 8) {print(gene[i])}
  if (!(tmp %in% count)) {
    count <- c(count, tmp)
  }
}
count <- c(0, count)
count <- sort(count)

#compute number of genes in each class ####
df <- data.frame(num_switch = count, num_gene = c(rep(0, length(count))))
df[df$num_switch == 0, ]$num_gene = length(unique(ref$gene_id)) - length(gene)
for (i in 1:length(gene)) {
  tmp <- sum(data$gene_id.y == gene[i])
  df[df$num_switch == tmp, ]$num_gene = df[df$num_switch == tmp, ]$num_gene + 1
}


#make barplot with gapped y axis ####
#whole figure
ggplot(data=df,aes(as.factor(num_switch),num_gene,fill=as.factor(num_switch)))+
  geom_bar(stat = "identity")+
  theme_classic()
  

#first fig
p1=ggplot(data=df,aes(as.factor(num_switch),num_gene,fill=as.factor(num_switch)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=num_gene), vjust=0, size = 5) + 
  coord_cartesian(ylim = c(0,800))+ #setting the down part
  theme_classic() +
  theme(legend.position="none",axis.title=element_text(size=16)) +
  xlab("Number of changes in transcript usage") +
  # ylab("Number of gene")
  ylab("Number of differentially expressed gene")
  

#second fig
p2=ggplot(data=df,aes(as.factor(num_switch),num_gene,fill=as.factor(num_switch)))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=num_gene), vjust=0, size = 5) + 
  coord_cartesian(ylim = c(12000,13000))+ #setting the upper part
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
pdf(paste0(plot_dir_out, 'isoswitch_per_gene.pdf'))
# pdf(paste0(plot_dir_out, 'isoswitch_per_DE_gene.pdf'), 8, 6)
p2/p1
dev.off()


#compare porportions of genes with transcript usage changes in all genes vs DE genes ####
df <- data.frame(DE = c(rep("Background", 2), rep("Differentially expressed", 2), rep("Not differentially expressed", 2)),
                 isoChange = c(rep(c("At least one", "No"), 3)), 
                 count = c(9291, length(unique(background_ref$gene_id))-9291,
                           6966, length(unique(DE_ref$gene_id))-6966,
                           2325, length(unique(NotDE_ref$gene_id))-2325))

df_fisherstest <- t(matrix(c(9291, length(unique(background_ref$gene_id))-9291,
         6966, length(unique(DE_ref$gene_id))-6966,
         2325, length(unique(NotDE_ref$gene_id))-2325), 2, 3))

row.names(df_fisherstest) <- c("Background", "Differentially expressed", "Not differentially expressed")
colnames(df_fisherstest) <- c("At least one", "No")


df_fisherstest <- as.table(df_fisherstest) %>%
  pairwise_fisher_test() %>%
  mutate(y.position = c(14800, 15600, 14800))
  # mutate(y.position = c(1.01, 1.06, 1.01))




pdf(paste0(plot_dir_out, 'bar_compare_genes_DTU_number.pdf'), 6, 4)
ggplot(df, aes(x = DE, y = count)) +
  geom_bar(aes(fill = isoChange), position = 'stack', stat = 'identity') +
  # geom_bar(aes(fill = isoChange), position = 'fill', stat = 'identity') +
  # scale_y_continuous(labels = scales::percent, limits = c(0,1.06,0.1), oob = rescale_none)+
  scale_x_discrete(labels = c("All", "Differentially\nexpressed", "Not differentially\nexpressed"))+
  # stat_pvalue_manual(df_fisherstest,label = "p.adj.signif", size = 3.7, tip.length = 0, bracket.shorten = 0.1,
  #                    coord.flip = TRUE) +
  theme_classic()+
  xlab("Gene types") +
  # ylab("Proportion of genes") +
  ylab("Number of genes") +
  # labs(title= "Pearson's Chi-squared test p-value < 0.001") + 
  guides(fill=guide_legend(title="Differential\ntranscript usage")) +
  theme(text = element_text(size = 14)) + 
  coord_flip()

dev.off()





###3.gene type of genes with isoswitch ####
tmp <- data %>% select(condition_1, condition_2, gene_name, gene_id.y, gene_type)
tmp <- tmp %>% mutate(transition = paste0(condition_1, '_vs_', condition_2))
pdf(paste0(plot_dir_out, 'bar_subsequent_gene_type.pdf'), 8, 6)
ggplot(tmp, aes(x = transition, fill = gene_type)) +
  geom_bar(position = 'stack', stat = 'count') +
  theme_bw()+
  xlab("Subsequent state transition") +
  ylab("Frequency of genes with isoform switch") +
  guides(fill=guide_legend(title="Gene type")) + 
dev.off()

###4.overlapping between isoswitch genes  ####
summary <- data.frame(group = rep("", 4), number = rep(0, 4))
#use unfiltered genes with significant isofrom switches between subsequent states ####
cond<-c()
genes<-list()
for (i in 1:(length(days)-1)){
  t<-paste0(days[i],' vs ',days[i+1])
  summary$group[i] <- t
  cond<-c(cond,t)
  # fn<-paste0(all_dire,'GenesWithIsoSwitch_day',days[i],'_vs_day',days[i+1],'.txt')
  fn<-paste0(all_dire,'GenesWithIsoSwitch_day',days[i],'_vs_day',days[i+1],'.csv')
  df<-read.csv(fn,header=TRUE,col.names = c('temp'))
  summary$number[i] <- nrow(df)
  genes<-append(genes,list(df$temp))
}
names(genes)<-cond

#output genes with isoform switch at any time
flat_genes <- unname(unlist(genes))
length(unique(flat_genes))
list <- as.data.frame(unique(flat_genes))
write.csv(list, file = paste0(res_dir_out, 'GenesWithIsoSwitch_at_any_time.csv'), row.names = FALSE, col.names = FALSE, quote = FALSE)


pdf(paste0(plot_dir_out,'Venn_changes_subsequent.pdf'),3,3)
venn(genes, ilab=TRUE, zcolor = "style", box = F, ilcs = 0.8, sncs = 1)
dev.off()

pdf(paste0(plot_dir_out,'gene_number_subsequent.pdf'),3,3)
ggplot(summary, aes(x = as.factor(group), y = number, fill = as.factor(group))) +
  geom_bar(stat = "identity") +
  # geom_text(aes(label = number), vjust = 0) + 
  scale_fill_hue(c = 40) +
  theme_classic() +
  theme(legend.position="none") +
  xlab("Comparison of consecutive days") +
  ylab("Number of genes with\n transcript isoform changes")
dev.off()

#check the common isoSwitch genes across all phase transitions
common <- Reduce(intersect, list(genes[1]$`0 vs 1`, genes[2]$`1 vs 2`, 
                       genes[3]$`2 vs 3`, genes[4]$`3 vs 5`))




#use filtered dominant-isoform-switch genes ####
cond<-c()
genes<-list()
for (i in 1:(length(days)-1)){
  t<-paste0(days[i],' vs ',days[i+1])
  summary$group[i] <- t
  cond<-c(cond,t)
  fn<-paste0(filt_dire,'mainIsoGenes_day',days[i],'_vs_day',days[i+1],'.txt')
  df<-read.csv(fn,header=FALSE,col.names = c('temp'))
  summary$number[i] <- nrow(df)
  genes<-append(genes,list(df$temp))
}
names(genes)<-cond

pdf(paste0(dirout,'filtered_Venn_changes_subsequent.pdf'),3,3)
venn(genes, ilab=TRUE, zcolor = "style", box = F)
dev.off()

pdf(paste0(dirout,'filtered_gene_number_subsequent.pdf'),3,4)
ggplot(summary, aes(x = as.factor(group), y = number, fill = as.factor(group))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = number), vjust = 0) + 
  scale_fill_hue(c = 40) +
  theme_bw()+
  theme(legend.position="none") +
  xlab("Phase") +
  ylab("Number of genes with isoform switch after filtering")
dev.off()


#use filtered NotDE & dominant-isoform-switch genes ####
cond<-c()
genes<-list()
for (i in 1:(length(days)-1)){
  t<-paste0(days[i],' vs ',days[i+1])
  summary$group[i] <- t
  cond<-c(cond,t)
  fn<-paste0(filt_dire,'mainIsoGenesNotDE_day',days[i],'_vs_day',days[i+1],'.txt')
  df<-read.csv(fn,header=FALSE,col.names = c('temp'))
  summary$number[i] <- nrow(df)
  genes<-append(genes,list(df$temp))
}
names(genes)<-cond

pdf(paste0(dirout,'filteredNotDE_Venn_changes_subsequent.pdf'),3,3)
venn(genes, ilab=TRUE, zcolor = "style", box = F)
dev.off()

pdf(paste0(dirout,'filteredNotDE_gene_number_subsequent.pdf'),3,4.5)
ggplot(summary, aes(x = as.factor(group), y = number, fill = as.factor(group))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = number), vjust = 0) + 
  scale_fill_hue(c = 40) +
  theme_bw()+
  theme(legend.position="none") +
  xlab("Phase") +
  ylab("Number of genes with isoform switch after filtering_NotDE")
dev.off()


###5. comparison of types of AS events during different transitions ####
event_df <- read.csv(event_fn, row.names = 1)
event_df <- event_df[event_df$Comparison == "day0 vs day1" | 
                       event_df$Comparison == "day1 vs day2" |
                       event_df$Comparison == "day2 vs day3" |
                       event_df$Comparison == "day3 vs day5" , ]
pdf(paste0(plot_dir_out, 'bar_subsequent_AS_type.pdf'), 8, 6)
ggplot(event_df, aes(x = Comparison, y = nrGenesWithConsequences, fill = AStype)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw()+
  xlab("Subsequent state transition") +
  ylab("Proportion of genes with isoform switch") +
  guides(fill=guide_legend(title="AS type"))
dev.off()

###6. genes with isoform switch in Volker's gene list
gene_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/GenesWithIsoSwitch_at_any_time.csv'
gene_df <- read.csv(gene_fn, header = 1)
colnames(gene_df)[1] <- 'gene_id'
gene_df <- gene_df %>% left_join(ref, by = 'gene_id')
gene_set <- gene_df$gene_name

goi_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/plotExpression/List/Volker_trafficking.txt'
goi_df <- read.csv(goi_fn, header = FALSE)
goi_set <- goi_df$V1

goi_switch <- intersect(gene_set, goi_set)
goi_switch <- as.data.frame(goi_switch)
write.csv(goi_switch, '/project/Neurodifferentiation_System/owlmayerTemporary/derek/plotExpression/List/Volker_trafficking_switching_RNAseq.txt', 
          row.names = FALSE, col.names = FALSE, quote = FALSE)


###6. Heatmap of cell cycle signatures from MSigDB (KEGG_CELL_CYCLE) ####
#loading the RSEM results
count_df <- read.csv("/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_gene_TPM_shortreads_byDay.csv", header = T)
day <- c('D0', 'D1', 'D2', 'D3', 'D5')
for (c in day){
  count_df[,paste0("Mean_",c)]<-count_df%>% select(starts_with(c))%>%rowMeans()
}

#load DE gene list ####
de <- read.csv("/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/DE_at_any_timepoint.txt", header = F)
for (i in 1:nrow(de)){
  de$V1[i] = unlist(strsplit(de$V1[i], '\\.'))[1]
}
de_count_df <- count_df  %>% filter(gene_id %in% de$V1)

#load cell cycle-related gene list ####
df <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/RawData/KEGG_CELL_CYCLE.v2023.1.Hs.grp", header = T)


ccTPM<- de_count_df  %>% filter(gene_symbol %in% df$KEGG_CELL_CYCLE)
rownames(ccTPM) <- ccTPM$gene_symbol
meanTPM <- ccTPM %>% select(starts_with("Mean"))
#z-scale the mean TPM for each transcript ####
standTPM<-as.data.frame(t(scale(t(meanTPM))))


fviz_nbclust(standTPM, kmeans, method = "wss")
fviz_nbclust(standTPM, kmeans, method = "silhouette")
fviz_nbclust(standTPM, kmeans, method = "gap_stat")

#kmeans with k = 2 ####
kclusters<-kmeans(standTPM,centers = 2, nstart=25)
standTPM$maxi<-max.col(standTPM)
standTPM$max<- sapply(as.data.frame(t(standTPM)),max)
standTPM$class<-kclusters$cluster
standTPM<- standTPM %>% arrange(desc(class))




kclasses<-standTPM %>%select(class)
kclasses$class <- factor(kclasses$class, levels = c(1,2))
standTPM<-standTPM%>%select(-c(class,max,maxi))

Var1        <- qualitative_hcl(2,palette="Dark3")
names(Var1) <- factor(c(1,2))
anno_colors <- list('class' = Var1)

pheatmap(standTPM,show_rownames=T,cluster_cols=F,cluster_rows = F, annotation_colors=anno_colors,
         annotation_row=kclasses, annotation_names_row = F,
         labels_col = c('day 0', 'day 1', 'day 2', 'day 3', 'day 5'), angle_col = '0',
         color =rev(scico(100,palette='roma'))[10:90],
         annotation_legend = F, fontsize_col = 12, fontsize_row = 6, border_color = NA, fontsize = 12, 
         gaps_row = head(as.numeric(cumsum(rev(table(kclasses$class)))), -1),
         filename = paste0(plot_dir_out, 'heatmap_cell_cycle.pdf'), width = 6, height = 8
)

###6. Heatmap of neuronal type markers ####
#load neuronal markers ####
marker <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/RawData/neural_marker.csv", header = T)
#focus on ####
test <- c("Immature neurons", "Mature neurons", "GABAergic neurons", "Glutamatergic neurons",
          "Dopaminergic neurons", "Serotonergic neurons", "Cholinergic neurons")
marker <- marker %>% filter(type %in% test)


marker$type <- factor(marker$type, test)

marker <- marker%>% arrange(factor(type, levels = test))
#load TPM values ####
count_df <- read.csv("/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_gene_TPM_shortreads_byDay.csv", header = T)
day <- c('D0', 'D1', 'D2', 'D3', 'D5')
for (c in day){
  count_df[,paste0("Mean_",c)]<-count_df%>% select(starts_with(c))%>%rowMeans()
}


marker_count <- marker %>% left_join(count_df, by = c("gene" = "gene_symbol"))
meanTPM <- marker_count %>% select(starts_with("Mean"))

## delete ADGRE1 from Microglia since TPM = 0 for all 5 days ####
# meanTPM <- meanTPM[-48, ]
# marker <- marker[-48, ]
#z-scale the mean TPM for each transcript ####
standTPM<-as.data.frame(t(scale(t(meanTPM))))
rownames(standTPM) <- paste0("row_", seq(nrow(standTPM)))


annot_row <- marker %>% select(type)
rownames(annot_row) <- paste0("row_", seq(nrow(annot_row)))



# Var1        <- qualitative_hcl(15,palette="Dark3")
Var1        <- qualitative_hcl(7,palette="Dynamic")
names(Var1) <- factor(unique(marker$type))
anno_colors <- list('type' = Var1)

pheatmap(standTPM,show_rownames=T,cluster_cols=F,cluster_rows = F, annotation_colors=anno_colors,
         annotation_row=annot_row, annotation_names_row = F,
         labels_row = marker$gene, 
         labels_col = c('day 0', 'day 1', 'day 2', 'day 3', 'day 5'), angle_col = '0',
         color =rev(scico(100,palette='roma'))[10:90],
         annotation_legend = T, fontsize_col = 12, fontsize_row = 8, border_color = NA, fontsize = 12, 
         gaps_row = head(as.numeric(cumsum(table(annot_row$type))), -1),
         filename = paste0(plot_dir_out, 'heatmap_neuron_type.pdf'), width = 8, height = 5
)


