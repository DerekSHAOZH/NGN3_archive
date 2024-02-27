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
library(data.table)
library(ggplotify)
library(reshape2)
library(ggnewscale)
library(scales)



plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Plots/'

###Example: rMAPs pval tables of day3 vs day2 exon skipping ####
df_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitches_RBPs/rMAPS/Results/day3_background_day2/SE/pVal.up.vs.bg.RNAmap.txt"
df <- read.table(df_fn, header = 1)

HuR_df <- df %>% filter(RBP %like% 'HuR')

###Heatmap of rMAPs pval ####
day <- c(0,1,2,3,5)
threshold=3

##only UP  ####
df_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/up.vs.bg.min_pval.csv"
df <- read.csv(df_fn)
rownames(df) <- df$RBP
data <- df[, -1]
data <- -log10(data)




# fviz_nbclust(data, kmeans, method = "wss")
# fviz_nbclust(data, kmeans, method = "silhouette")
# fviz_nbclust(data, kmeans, method = "gap_stat")
# 
# kclusters<-kmeans(data,centers = 4, nstart=25)
# data$class<-kclusters$cluster
# data<- data %>% arrange(desc(class))
# kclasses<-data %>%select(class)

kclasses <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/RBP.SE_day3_vs_day2.cluster_result.csv",
                     row.names = 1)
kclasses$class <- factor(kclasses$class, levels = c(1,2,3,4))

# write.csv(kclasses, "/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/RBP.SE_day3_vs_day2.cluster_result.csv",
#           quote = F)


# data<-data%>%select(-c(class))
data <- as.matrix(data)
ind <- which(data <= threshold)
data[ind] <- NA

data <- as.data.frame(data)
data <- data %>% filter(!is.na(day1.vs.day0) | !is.na(day2.vs.day1) | !is.na(day3.vs.day2) | !is.na(day5.vs.day3))
sig <- rownames(data)
data <- as.matrix(data)

kclasses <- as.data.frame(kclasses[sig, ])
rownames(kclasses) <- sig
colnames(kclasses) <- "group"
kclasses <- kclasses[order(kclasses$group, decreasing = T), , drop  =F]
data <- data[rownames(kclasses), ]
# colnames(data) <- c("day1 vs day0","day2 vs day1", "day3 vs day2",  "day5 vs day3")

#heatmap: time-specific ####
breaksList = seq(3, 10, by = 0.01)
cols <- c("#0079FF", "#00DFA2",  "#FF7F3F", "#FF0060")
names(cols) <- factor(c(1,2,3,4))
anno_colors <- list('group' = cols)
# pheatmap(data,
#          # show_rownames=T,cluster_rows=F, 
#          show_rownames=F,cluster_rows=F, 
#          # show_colnames=T,cluster_cols=F,
#          show_colnames=F,cluster_cols=F,
#          annotation_colors=anno_colors,
#          annotation_row=kclasses, annotation_names_row = F,
#          # annotation_col=annot_col, annotation_names_col = F, 
#          labels_col = c(expression("day 0" %->% "day 1"), expression("day 1" %->% "day 2"), 
#                         expression("day 2" %->% "day 3"), expression("day 3" %->% "day 5")), 
#          angle_col = '45',
#          color = colorRampPalette( brewer.pal(9, "Reds") )(length(breaksList)),
#          annotation_legend = F,
#          fontsize = 14, border_color = NA, 
#          breaks = breaksList, 
#          gaps_row = head(as.numeric(cumsum(rev(table(kclasses$group)))), -1),
#          # filename = paste0(plot_dir_out, 'heatmap_rmaps_pval.pdf'), width = 8, height = 12
# )

melt_data <- reshape2::melt(as.matrix(data))

annot <- kclasses %>% mutate(rbp = rownames(kclasses), x = "1")
melt_data <- melt_data %>% mutate(group = rep(annot$group, 4))



melt_data$Var1 <- factor(melt_data$Var1, levels=rev(rownames(data)))
melt_data$group <- factor(melt_data$group, levels=rev(c(1,2,3,4)))
annot$rbp <- factor(annot$rbp, levels=rev(annot$rbp))
annot$group <- factor(annot$group, levels=rev(c(1,2,3,4)))

cooccur_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/cooccurrence_pvalue_motif_clip.txt"
cooccur_df <- read.delim(cooccur_fn)

cooccur_df$comparison <- sub("day(\\d+)_background_day(\\d+)", "day\\1.vs.day\\2", cooccur_df$comparison)
cooccur_df <- cooccur_df %>% filter(diff == "up")
cooccur_df <- cooccur_df %>% left_join(annot, by = "rbp")


meta_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Data/RBP_peak_meta.txt"
meta_df <- read.table(meta_fn, sep = '\t', header = F)
annot <- annot %>% left_join(meta_df[, c("V1", "V3")], by = c("rbp" = "V1"))
annot[is.na(annot$V3), "V3"] <- "No"
annot$V3 <-  factor(annot$V3, levels=c("Yes", "No"))
annot$x_1 <- "2"


breaksList = seq(3, 10, by = 0.01)
reds <- colorRampPalette( brewer.pal(9, "Reds") )(length(breaksList))
sig_heatmap <- ggplot(melt_data, aes(x = Var2, y = Var1)) +
  # Important for ggnewscale is to specify a fill in the layer/geom itself
  geom_tile(aes(fill = value),
             colour = "transparent") +
  # geom_point(data = df[df$pval == "FDR<0.05",]) +
  scale_fill_gradientn(colours = reds,
                       limits=c(3,10),
                       na.value="#DDDDDD",
                       name = "Consensus motif enrichment",
                       oob=scales::squish) +
  # Set new scale fill after you've specified the scale for the heatmap
  new_scale_fill() +
  geom_tile(data = annot, aes(x, rbp, fill = group),
            width = 0.3, colour = "transparent") +
  scale_fill_manual(name = "Class", values = cols) +
  new_scale_fill() +
  geom_tile(data = annot, aes(x_1, rbp, fill = V3),
            width = 0.3, colour = "transparent") +
  scale_fill_manual(name = "CLIP-seq data", values = c("Yes" = "#492E87",
                                                       "No" = "#F3CCF3")) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y")+
  # scale_x_discrete(name = "", expand = c(0,0)) +
  # scale_y_discrete(name = "", expand = c(0,0),
  #                  limits = c(levels(df$Var1), "Condition"),
  #                  position = "right") +
  # coord_equal() +
  geom_point(data = cooccur_df, mapping = aes(x = comparison, y = rbp, color = "p < 0.05"),
             shape = 23, fill = "white", size = 1, stroke = 1) +
  scale_color_manual(name = "CLIP-seq binding site enrichment",
                     values = c("p < 0.05" = "black")) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        # axis.text.x = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "left") 
  



#rbp name panel ####
name_df <- data.frame("rbp" = rownames(data),
                      "group" = kclasses$group,
                   "pos" = 1)
name_df$rbp <- factor(name_df$rbp, levels=rev(name_df$rbp))
name_df$group <- factor(name_df$group, levels=rev(c(1,2,3,4)))
bg_df <- name_df %>%
  mutate(color = ifelse(row_number() %% 2 == 1, "lightgrey", "transparent"))
name_panel <- ggplot() +
  geom_tile(data = bg_df, aes(y = rbp, x = pos, fill = color), 
            alpha = 0.5, inherit.aes = FALSE) +
  geom_text(name_df, mapping = aes(x = pos, y = rbp, label = rbp), size = 2.5) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") + 
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank()
        
  ) +
  scale_fill_identity()  

#rbp gene expression ####
exp_df <- data.frame("rbp" = rownames(data),
                     "group" = kclasses$group)
exp_df$rbp <- factor(exp_df$rbp, levels=rev(exp_df$rbp))
exp_df$group <- factor(exp_df$group, levels=rev(c(1,2,3,4)))

meta_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Data/RBP_peak_meta.txt"
meta_df <- read.table(meta_fn, sep = '\t', header = F)

exp_df <- exp_df %>% left_join(meta_df, by = c("rbp" = "V1"))

#Use normalized counts as expression values
info_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_gene_TPM_shortreads_byDay.csv"
info_df <- read.csv(info_fn)
exp_df <- exp_df %>% left_join(info_df[, c("gene_id", "gene_symbol")], by = c("V2" = "gene_symbol"))

gene_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/counts_normalized.txt"
gene_df <- read.delim(gene_fn, header = T)
gene_df$gene_id <- sub("\\..*$", "", rownames(gene_df))
exp_df <- exp_df %>% left_join(gene_df, by = "gene_id")

#Use TPM as expression values
# gene_fn <- "/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_gene_TPM_shortreads_byDay.csv"
# gene_df <- read.csv(gene_fn)
# 
# colnames(gene_df)[2:25] <- sub("D(\\d+)_(\\d+)", "day\\1_\\2", colnames(gene_df)[2:25])
# exp_df <- exp_df %>% left_join(gene_df, by = c("V2" = "gene_symbol"))


exp_df <- exp_df %>% select(contains("day"))
exp_df <- log10(exp_df+1)
# exp_df <- log2(exp_df+1)
# exp_df<-as.data.frame(t(scale(t(exp_df))))
rownames(exp_df) <- rownames(data)
condition <- sub("_.*$", "", colnames(exp_df))

# exp_heatmap <- pheatmap(exp_df,
#          show_rownames=F,cluster_rows = F,
#          show_colnames = F, cluster_cols=F, 
#          border_color = "transparent", 
#          color =rev(scico(100,palette='roma'))[10:90],
#          gaps_row = head(as.numeric(cumsum(rev(table(kclasses$group)))), -1),
#          gaps_col = head(as.numeric(cumsum(table(condition))), -1)
#          )
# exp_heatmap <- as.ggplot(exp_heatmap)

melt_exp_df <- reshape2::melt(as.matrix(exp_df))

melt_exp_df <- melt_exp_df %>% left_join(name_df, by = c("Var1" = "rbp"))
melt_exp_df$Var1 <- factor(melt_exp_df$Var1, levels=rev(name_df$rbp))
melt_exp_df$group <- factor(melt_exp_df$group, levels=rev(c(1,2,3,4)))

melt_exp_df <- melt_exp_df %>% mutate(condition = sub("_.*$", "", Var2))
melt_exp_df$condition <- factor(melt_exp_df$condition, levels=c("day0", "day1", "day2", "day3", "day5"))

break1 <- 

exp_heatmap <- ggplot(melt_exp_df, aes(x = Var2, y = Var1)) +
  # Important for ggnewscale is to specify a fill in the layer/geom itself
  geom_tile(aes(fill = value),
            colour = "transparent") +
  # geom_point(data = df[df$pval == "FDR<0.05",]) +
  # scale_fill_gradientn(colours = rev(scico(100,palette='roma'))[10:90],
  #                      # limits=c(-2,2),
  #                      na.value="#DDDDDD",
  #                      # name = "Z-normalized RLE",
  #                      name = "log10(RLE+1)",
  #                      oob=scales::squish) +
  scale_fill_gradientn(colours = c(rev(scico(136,palette='roma'))[1:68], rev(scico(54,palette='roma'))[28:54]),
                       # values = rescale(c(min(as.matrix(exp_df), na.rm = T), mean(as.matrix(exp_df), na.rm = T), max(as.matrix(exp_df), na.rm = T))),
                       values = rescale(c(seq(min(as.matrix(exp_df), na.rm = T), median(as.matrix(exp_df), na.rm = T), by = 0.05), seq(median(as.matrix(exp_df), na.rm = T), max(as.matrix(exp_df), na.rm = T), by = 0.05))),
                       # limits=c(-2,2),
                       na.value="#DDDDDD",
                       # name = "Z-normalized RLE",
                       name = "log10(RLE+1)",
                       oob=scales::squish) +
  facet_grid(group ~ condition, scales = "free", space = "free")+
  # scale_x_discrete(name = "", expand = c(0,0)) +
  # scale_y_discrete(name = "", expand = c(0,0),
  #                  limits = c(levels(df$Var1), "Condition"),
  #                  position = "right") +
  # coord_equal() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        # axis.text.x = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "right") 



# pdf(file = paste0(plot_dir_out, "combined_heatmap.log10RLE.pdf"), 9,6)
ggarrange(sig_heatmap, name_panel, exp_heatmap, nrow = 1, ncol = 3,
          widths = c(1.05,0.15,0.6))
dev.off()




##UP & DN ####

up_df_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/up.vs.bg.min_pval.csv"
up_df <- read.csv(up_df_fn)
rownames(up_df) <- up_df$RBP
up_data <- up_df[, -1]
up_data <- -log10(up_data)

# rownames(up_data) <- paste0(rownames(up_data), "_up")
colnames(up_data) <- paste0(colnames(up_data), "_up")

dn_df_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/dn.vs.bg.min_pval.csv"
dn_df <- read.csv(dn_df_fn)
rownames(dn_df) <- dn_df$RBP
dn_data <- dn_df[, -1]
dn_data <- log10(dn_data)

# rownames(dn_data) <- paste0(rownames(dn_data), "_dn")
colnames(dn_data) <- paste0(colnames(dn_data), "_dn")

# data <- rbind(up_data, dn_data)
data <- cbind(up_data, dn_data)

# data <- data[order(rownames(data)), ]
data <- data[, order(colnames(data))]

# fviz_nbclust(data, kmeans, method = "wss")
# fviz_nbclust(data, kmeans, method = "silhouette")
# fviz_nbclust(data, kmeans, method = "gap_stat")
# 
# kclusters<-kmeans(data,centers = 3, nstart=25)
# data$class<-kclusters$cluster
# data<- data %>% arrange(desc(class))
# kclasses<-data %>%select(class)

kclasses <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/RBP.SE_day3_vs_day2.cluster_result.csv",
                     row.names = 1)
# kclasses$class <- factor(kclasses$class, levels = c(1,2,3, 4))
kclasses$class <- factor(kclasses$class, levels = c(1,2,3))

# write.csv(kclasses, "/project/Neurodifferentiation_System/owlmayerTemporary/derek/motif_rMAPs/Results/RBP.SE.cluster_result.up_and_dn.csv",
#           quote = F)


# data<-data%>%select(-c(class))


data <- as.matrix(data)
ind <- which(abs(data) <= threshold)
data[ind] <- NA

data <- as.data.frame(data)
# data <- data %>% filter(!is.na(day1.vs.day0) | !is.na(day2.vs.day1) | !is.na(day3.vs.day2) | !is.na(day5.vs.day3))
data <- data %>% filter(if_any(everything(), ~ !is.na(.)))
data <- as.matrix(data)

sig <- rownames(data)

kclasses <- as.data.frame(kclasses[sig, ])
rownames(kclasses) <- sig
colnames(kclasses) <- "group"

# annot_row <- data.frame(do.call(rbind, strsplit(rownames(data), "_")))
# colnames(annot_row) <- c("RBP", "Differential splicing")
# rownames(annot_row) <- rownames(data)
# annot_row[annot_row[,2] == "up", 2] <- "Included"
# annot_row[annot_row[,2] == "dn", 2] <- "Skipped"
# annot_row_new <- annot_row[, -1, drop = F]

annot_col <- data.frame(do.call(rbind, strsplit(colnames(data), "_")))
colnames(annot_col) <- c("comparison", "Differential splicing")
rownames(annot_col) <- colnames(data)
annot_col[annot_col[,2] == "up", 2] <- "Included"
annot_col[annot_col[,2] == "dn", 2] <- "Skipped"
annot_col_new <- annot_col[, -1, drop = F]


#heatmap: time-specific ####
up_breaksList = seq(3, 10, by = 0.01)
dn_breaksList = seq(-10, -3, by = 0.01)
colors <- c(colorRampPalette( rev(brewer.pal(9, "Blues")) )(length(dn_breaksList)),
            colorRampPalette( brewer.pal(9, "Reds") )(length(up_breaksList)))

cols <- c("#F8766D", "#619CFF")
names(cols) <- factor(c("Included","Skipped"))
# class_cols <- c("#0079FF", "#00DFA2",  "#FF7F3F", "#FF0060")
class_cols <- c("#0079FF", "#00DFA2",  "#FF7F3F")
# names(class_cols) <- factor(c(1,2,3,4))
names(class_cols) <- factor(c(1,2,3))


anno_colors <- list("Differential splicing" = cols, 
                    'group' = class_cols)
pheatmap(data,
         show_rownames=T,cluster_rows=F, 
         show_colnames=F,cluster_cols=F,
         annotation_colors=anno_colors,
         annotation_row=kclasses, annotation_names_row = F,
         annotation_col=annot_col_new, annotation_names_col = F,
         labels_col = c(expression("day 0" %->% "day 1"), expression("day 1" %->% "day 2"), 
                        expression("day 2" %->% "day 3"), expression("day 3" %->% "day 5")), 
         angle_col = '45',
         color = colors,
         # annotation_legend = T, 
         fontsize = 14, border_color = NA, 
         breaks = c(dn_breaksList, up_breaksList),
         gaps_row = head(as.numeric(cumsum(rev(table(kclasses$group)))), -1),
         gaps_col = head(as.numeric(cumsum(rev(table(annot_col$comparison)))), -1),
         filename = paste0(plot_dir_out, 'heatmap_rmaps_pval.without_col_label.pdf'), width = 8, height = 12
)

dev.off()


#test ####
data <- read.delim(
  system.file("demo-data/housetasks.txt", package = "ggpubr"),
  row.names = 1
)
data

# test <- data.frame("x_x" = c(4,3), 
#                    "y_y" = c(1,2))
test <- data.frame("x_x" = "Wife", 
                   "y_y" = "Laundry")
ggballoonplot(data, ggtheme = theme_minimal(base_size = 14)) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.title.position = "plot",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.ticks = element_blank()) +
  geom_tile(color = "black", fill = "transparent") +
  geom_point(data = test, mapping = aes(x = x_x, y = y_y)) +
  coord_equal()
 

  

