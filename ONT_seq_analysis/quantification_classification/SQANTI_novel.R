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
library(scales) # to access break formatting functions
plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Plots/'
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)


##load SQANTI classification result ####
# class_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_classification_TPM.txt'
# class_df <- read.csv(class_fn, sep = '\t')
#after filtering
class_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter.1/nanopore_RulesFilter_result_classification.txt'
# class_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_RulesFilter_result_classification.txt'
class_df <- read.csv(class_fn, sep = '\t')
class_df <- class_df %>% filter(filter_result == "Isoform")

# info_df <- class_df %>% select("isoform", "chrom", "strand", starts_with("associated"))
# write.table(info_df, "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore.filtered.transcript_novelty.txt",
            # quote = F, row.names = F, sep = '\t')
##load SQANTI junction ####
junction_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_QC/nanopore_junctions.txt'
junction_df <- read.csv(junction_fn, sep = '\t')


###Distribution of active novel transcripts (NIC, NNC, fusion) per day ####
###"active" criteria: filter out transcripts with mean TPM < 5, then mean TPM across replicates per day > 5
class_df <- class_df %>% rowwise() %>% mutate(FL.mean = mean(c(FL.day0_1, FL.day0_2, FL.day0_3, FL.day3_1, FL.day3_2, FL.day3_3,
                                       FL.day5_1, FL.day5_2, FL.day5_3), na.rm = T),
                    FL.day0_mean = mean(c(FL.day0_1, FL.day0_2, FL.day0_3), na.rm = T),
                    FL.day3_mean = mean(c(FL.day3_1, FL.day3_2, FL.day3_3), na.rm = T),
                    FL.day5_mean = mean(c(FL.day5_1, FL.day5_2, FL.day5_3), na.rm = T))

active <- class_df[class_df$FL.mean > 5 & !is.na(class_df$FL.mean), ]

df = data.frame(day = c(rep('day 0', 5), rep('day 3', 5), rep('day 5',5)),
                type = c(rep(c('FSM', 'ISM', 'NIC', 'NNC', 'Fusion'), 3)),
                num = rep(0,15))
df$type <- factor(df$type, levels = c("FSM", "ISM", 
                                        "NIC", "NNC", 'Fusion'))

day <- c('day0', 'day3', 'day5')
type <- c('full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'fusion')

for (i in 1:length(day)){
  for (j in 1:length(type)){
    day_ind <- which(colnames(active) == paste0('FL.', day[i], '_mean'))
    tmp <- sum((active$structural_category == type[j]) & (active[, day_ind]>5))
    df$num[(i-1)*5+j] <- tmp
  }
}

cols = c("FSM" = "#6baed6", "ISM" = "#0EA293", 
         "NIC" = "#fc8d59", "NNC" = "#ee6a50",
         'Fusion' = '#FF2171')

#all types: proportion ####
pdf(paste0(plot_dir_out, 'long_active_transcript_type_proportion_per_day.pdf'), 8, 4)
ggplot(df, aes(x = day, y= num, fill = type)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw()+
  xlab("Day") +
  ylab("Proportion of active transcripts") +
  guides(fill=guide_legend(title="Transcript isoform types")) +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=16))
  # theme(axis.text = element_text(size = 20))
dev.off()

#all types: frequency ####
pdf(paste0(plot_dir_out, 'long_active_transcript_type_frequency_per_day.pdf'), 8, 4)
ggplot(df, aes(x = day, y= num, fill = type)) +
  geom_bar(position = 'stack', stat = 'identity') +
  theme_bw()+
  xlab("Day") +
  ylab("Number of active transcripts") +
  guides(fill=guide_legend(title="Transcript isoform types")) +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=16))
dev.off()




#only novel: proportion ####
pdf(paste0(plot_dir_out, 'novel_active_transcript_type_proportion_per_day.pdf'), 8, 4)
ggplot(df[df$type == 'NIC' | df$type == 'NNC' | df$type == 'Fusion', ], aes(x = day, y= num, fill = type)) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw()+
  xlab("Day") +
  ylab("Proportion of active transcript isoforms") +
  guides(fill=guide_legend(title="Transcript isoform types")) +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=16))
dev.off()

#only novel: frequency ####
pdf(paste0(plot_dir_out, 'novel_active_transcript_type_frequency_per_day.pdf'), 8, 4)
ggplot(df[df$type == 'NIC' | df$type == 'NNC' | df$type == 'Fusion', ], aes(x = day, y= num, fill = type)) +
  geom_bar(position = 'stack', stat = 'identity') +
  theme_bw()+
  xlab("Day") +
  ylab("Number of active transcript isoforms") +
  guides(fill=guide_legend(title="Transcript isoform types")) +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=16))
dev.off()






###AS mechanisms of novel isoform formation ####
#load AS detail of ONT transcripts (onlysSitchingGenes = True) ####
as_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/AS_type_detail.csv'
as_df <- read.csv(as_fn)

novel <- class_df[class_df$structural_category == 'novel_in_catalog' | 
                    class_df$structural_category == 'novel_not_in_catalog' |
                    class_df$structural_category == 'fusion', ]
novel <- novel[, c(1, 6, 7)]
as_df <- as_df %>% left_join(novel, by = c('isoform_id' = 'isoform'))
as_df <- as_df[!is.na(as_df$structural_category), ]
## Novel transcripts in switching genes: 2934 ####
as_type <- c('ES', 'MEE', 'MES', 'IR', 'A5', 'A3', 'ATSS', 'ATTS')
df = data.frame(AStype = as_type, 
                num = 0)

for (i in 1:length(as_type)){
  df[i, 2] <- sum(as_df[, which(colnames(as_df) == as_type[i])])
}

df <- df %>% 
  mutate(perc = `num` / sum(`num`)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc, accuracy = 0.1))


df$display_value <- sapply(df$num, function(x){
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

cols = c("ES" = "#9467bdff", "MEE" = "#e377c2ff", 
         "MES" = "#7f7f7fff", "IR" = "#8c564bff",
         'A5' = '#ff7f0eff', 'A3' = '#1f77b4ff',
         'ATSS' = '#2ca02cff', 'ATTS' = '#d62728ff') 
#pie chart ####
pdf(paste0(plot_dir_out, 'novel_AS_mechanism_sqanti3_legend.pdf'), 8, 8)
ggplot(df, aes(x = "", y = num, fill = AStype)) +
  geom_col(color = "white") +
  geom_text(aes(label = paste0(labels, '\n(', display_value, ')')), color = c("white"),
            position = position_stack(vjust = 0.5),
            show.legend = FALSE,size = 7
  ) +
  geom_text(aes(1.6, label = AStype), color = c("black"),
            position = position_stack(vjust = 0.5), hjust = "outward",
            show.legend = FALSE,size = 7
  ) +
  guides(fill = guide_legend(title = "Alternative splicing mechanism", size = 10)) +
  scale_fill_manual(values = cols) +
  coord_polar(theta = "y", clip = "off") + 
  theme_void() +
  theme(legend.position = "none")
  
  # theme(legend.title=element_text(size=20),
  #       legend.text=element_text(size=20))
dev.off()

###Example of novel isoforms ####
exp <- active[order(active$FL.mean, decreasing = T), ]
type <- c('incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'fusion')
exp <- exp %>%filter(structural_category %in% type)
exp <- exp[exp$min_cov > 10 & exp$perc_A_downstream_TTS < 60 & 
             exp$within_CAGE_peak == TRUE & exp$within_polyA_site == TRUE & 
             exp$polyA_motif_found == TRUE &
             exp$diff_to_gene_TSS > -50 & exp$diff_to_gene_TSS < 50 &
             exp$diff_to_gene_TTS > -50 & exp$diff_to_gene_TTS < 50 &
             exp$all_canonical == "Canonical", ]

conseq_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/consequence_detailed_sqanti.csv'
conseq <- read.csv(conseq_fn)
isoup <- conseq$isoformUpregulated
isodown <- conseq$isoformDownregulated
exp_conseq <- exp[exp$isoform %in% isoup | exp$isoform %in% isodown, ]

###Number of isoforms before & after filtering ####
class_fn <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore_RulesFilter_result_classification.txt'
class_df <- read.csv(class_fn, sep = '\t')


type <- c('full-splice_match', 'incomplete-splice_match', 'novel_in_catalog', 'novel_not_in_catalog', 'fusion')
type_new <- c("Full Splice Match", "Incomplete Splice Match", 
               "Novel in Catalog", "Novel Not in Catalog", 'Fusion')
class_df <- class_df[class_df$structural_category %in% type, ]
for (i in 1:length(type)) {
  class_df[class_df$structural_category == type[i], ]$structural_category <- type_new[i]
}
class_df$structural_category <- factor(class_df$structural_category, levels = type_new)
class_df$filter_result <- factor(class_df$filter_result, levels = rev(c("Isoform", "Artifact")))

cols = c('Artifact'= '#FF7D7D', 'Isoform'= '#59CE8F')
cols = c("Full Splice Match" = "#6baed6", "Incomplete Splice Match" = "#0EA293", 
         "Novel in Catalog" = "#fc8d59", "Novel Not in Catalog" = "#ee6a50",
         'Fusion' = '#FF2171')

pdf(paste0(plot_dir_out, 'count_long_transcript_type_after_filtering.ver2.pdf'), 8, 4)
# ggplot(class_df, aes(x = structural_category, fill = filter_result)) +
#   geom_bar(position = 'stack', stat = 'count') +
#   theme_bw()+
#   xlab("Structural category") +
#   ylab("Number of transcript isoforms") +
#   # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#   #               labels = trans_format("log10", math_format(10^.x))) +
#   # annotation_logticks(sides = 'b')  + 
#   guides(fill=guide_legend(title="Filter result")) +
#   scale_fill_manual(values = cols) +
#   theme(text=element_text(size=16)) +
#   coord_flip()
ggplot(class_df, aes(x = filter_result, fill = structural_category)) +
  geom_bar(position = 'dodge', stat = 'count') +
  theme_bw()+
  xlab("Filter result") +
  ylab("Number of transcript isoforms") +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  # annotation_logticks(sides = 'b')  +
  guides(fill=guide_legend(title="Structural category", nrow = 2)) +
  scale_fill_manual(values = cols, 
                    limits = rev(levels(class_df$structural_category))) +

  theme(text=element_text(size=16), 
        legend.position = 'top',
        axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  coord_flip()
dev.off()

pdf(paste0(plot_dir_out, 'count.log10_long_transcript_type_after_filtering.pdf'), 8, 4)
ggplot(class_df, aes(x = structural_category, fill = filter_result)) +
  geom_bar(position= 'dodge', stat = 'count') +
  theme_bw()+
  xlab("Structural category") +
  ylab("Number of transcript isoforms") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = 'b') + 
  guides(fill=guide_legend(title="Filter result")) +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=16)) +
  coord_flip()
dev.off()



pdf(paste0(plot_dir_out, 'perc_long_transcript_type_after_filtering.pdf'), 8, 4)
ggplot(class_df, aes(x = structural_category, fill = filter_result)) +
  geom_bar(position = 'fill', stat = 'count') +
  theme_bw()+
  xlab("Structural category") +
  ylab("Percentage of transcript isoforms") +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  # annotation_logticks(sides = 'b')  + 
  guides(fill=guide_legend(title="Filter result")) +
  scale_fill_manual(values = cols) +
  theme(text=element_text(size=16)) +
  coord_flip()
dev.off()

###Number of novel isoforms peak at different days ####
type <- c('novel_in_catalog', 'novel_not_in_catalog', 'fusion')

#day 0 ####
day0_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/valid_transcript.peak_at_day0.csv"
day0_df <- read.csv(day0_fn, header = F)
colnames(day0_df) <- "transcript_id"
day0_df <- day0_df %>% left_join(class_df, by = c("transcript_id" = "isoform"))

novel_day0_df <- day0_df %>% filter(structural_category %in% type)
#Number of novel isoforms peak at day 0: 3434 ####
#day 3 ####
day3_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/valid_transcript.peak_at_day3.csv"
day3_df <- read.csv(day3_fn, header = F)
colnames(day3_df) <- "transcript_id"
day3_df <- day3_df %>% left_join(class_df, by = c("transcript_id" = "isoform"))

novel_day3_df <- day3_df %>% filter(structural_category %in% type)
#Number of novel isoforms peak at day 3: 3336 ####
#day 5 ####
day5_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/valid_transcript.peak_at_day5.csv"
day5_df <- read.csv(day5_fn, header = F)
colnames(day5_df) <- "transcript_id"
day5_df <- day5_df %>% left_join(class_df, by = c("transcript_id" = "isoform"))

novel_day5_df <- day5_df %>% filter(structural_category %in% type)
#Number of novel isoforms peak at day 5: 5244 ####


