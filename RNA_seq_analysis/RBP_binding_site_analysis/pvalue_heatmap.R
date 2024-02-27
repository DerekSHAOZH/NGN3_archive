library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(colorspace)
library(RColorBrewer)
library(ggrepel)
library(reshape2)



res_dir <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_RBP/"
plot_out_dir <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Plots/"

###create clip pvalue table ####
rbps <- list.files(res_dir)



ranges <- list(list(1:50), list(51:300), list(301:550), list(551:600), list(601:650), list(651:900), list(901:1150), list(1151:1200))
comparison <- c("day1_background_day0", "day2_background_day1", "day3_background_day2", "day5_background_day3")
diff <- c("up", "dn")

for (comp in comparison) {
  out_dir <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_pvalue_table/", comp)
  if(!dir.exists(out_dir)){
    dir.create(out_dir)
  }
  
  for (d in diff) {
    pvalue_df <- data.frame("RBP" = rbps,
                            "smallest_p_in_upstreamExon-3prime" = 0,
                            "smallest_p_in_upstreamExonIntron" = 0,
                            "smallest_p_in_upstreamIntron" = 0,
                            "smallest_p_in_targetExon-5prime" = 0,
                            "smallest_p_in_targetExon-3prime" = 0,
                            "smallest_p_in_downstreamIntron" = 0,
                            "smallest_p_in_downstreamExonIntron" = 0,
                            "smallest_p_in_downstreamExon-5prime" = 0)
    
    for (rbp in rbps) {
      ind <- paste0(rbp, "_", comp, ".", d)
      raw_fn <- paste0(res_dir, rbp, "/", comp, "/", ind, ".SE.MATS.JC.txt.pvalues.txt")
      raw_df <- read.table(raw_fn, sep = '\t', header = F)
      
      for (i in 1:length(ranges)){
        #in raw_df, p value is -log10 transformed
        smallest_p <- max(raw_df[ranges[[i]][[1]], 2], na.rm = T)
        pvalue_df[pvalue_df$RBP == rbp, i+1] <- smallest_p
      }
    }
    out_fn <- paste0(out_dir, "/", "RBP.pVal.", d,".vs.bg.txt")
    write.table(pvalue_df, out_fn, sep = '\t', quote= F, row.names = F)
  }
}

###motif vs clip pvalue analysis ####
threshold <- 0.001  
threshold_log <- -log10(0.05)

comparison <- c("day1_background_day0", "day2_background_day1", "day3_background_day2", "day5_background_day3")
diff <- c("up", "dn")

comp <- "day5_background_day3"
d <- "dn"
##pvalue for predicted motif enrichment ####
motif_fn <- paste0("/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitches_RBPs/rMAPS/Results/", comp, "/SE/pVal.", d, ".vs.bg.RNAmap.txt")
motif_df <- read.table(motif_fn, sep = '\t', header = T)

regions <- colnames(motif_df)[2:ncol(motif_df)]

# all_regions <- paste(rep(regions, each = 2), c("motif", "clip"), sep = "_")



motif_df <- motif_df %>%
  filter(rowSums(select(., contains('smallest_p_in_')) < threshold) > 0)


motif_df <- motif_df %>%
  mutate(RBP = sub("\\..*", "", RBP)) %>%
  arrange(RBP) %>%
  group_by(RBP) %>%
  summarise_all(.funs = list(min))

motif_df[1:nrow(motif_df), 2:ncol(motif_df)] <- -log10(motif_df[1:nrow(motif_df), 2:ncol(motif_df)])

# colnames(motif_df)[2:ncol(motif_df)] <- paste0(colnames(motif_df)[2:ncol(motif_df)], "_motif")
signif_rbp <- motif_df$RBP




##pvalue for clip binding site enrichment ####
clip_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_pvalue_table/", comp, "/RBP.pVal.", d, ".vs.bg.txt")
clip_df <- read.table(clip_fn, sep = '\t', header = T)

clip_df <- clip_df %>%
  filter(rowSums(select(., contains('smallest_p_in_')) > threshold_log) > 0)

meta_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Data/RBP_peak_meta.txt"
meta_df <- read.table(meta_fn, sep = '\t', header = F)

clip_df <- clip_df %>% left_join(meta_df, by = c("RBP" = "V2"))

rownames(clip_df) <- clip_df$V1
clip_data <- as.matrix(clip_df %>% select(contains("smallest")))

indices <- which(clip_data > -log10(0.05), arr.ind = TRUE)


clip_signif <- data.frame("pos" = colnames(clip_data)[indices[, 2]],
                          "rbp" = rownames(clip_data)[indices[, 1]])
clip_signif <- clip_signif %>% filter(rbp %in% signif_rbp)




# motif_df <- motif_df %>% filter(!is.na(V2))

# colnames(clip_df)[2:ncol(clip_df)] <- paste0(colnames(clip_df)[2:ncol(clip_df)], "_clip")
   

##with all motif RBP ####
# df <- motif_df %>% left_join(clip_df, by = c("V2" = "RBP"))
# 
# rownames(df) <- df$RBP
# data <- df[, all_regions]
# rownames(data) <- df$RBP

#correlation analysis ####
# motif_p <- df %>% select(contains('motif'))
# rownames(motif_p) <- paste0(df$RBP, "_motif")
# clip_p <- df %>% select(contains('clip'))
# rownames(clip_p) <- paste0(df$RBP, "_clip")
# 
# cor_mat <- cor(t(motif_p), t(clip_p), method = "spearman")
#   
#   
# cor_list <- diag(cor_mat)
# # which(cor_list == max(cor_list, na.rm = T))
# which(cor_list >0.5)
# df$RBP[which(cor_list >0.5)]

# cor_df <- data.frame("RBP" = df$RBP,
#                       "comparison" = comp,
#                       "diff" = d,
#                       "corr" = cor_list)
# out_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_pvalue_table/", comp, "/RBP.", d, ".pvalue.spearman_correlation.txt")
# 
# write.table(cor_df, out_fn, sep = '\t', quote = F, row.names = F)
p_cor_df <- data.frame()

for (comp in comparison) {
  for (d in diff) {
    tmp_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_pvalue_table/", comp, "/RBP.", d, ".pvalue.spearman_correlation.txt")
    tmp_df <- read.table(tmp_fn, sep = '\t', header = T)
    
    p_cor_df <- rbind(p_cor_df, tmp_df)
  }
}
p_cor_df$comparison <- as.factor(p_cor_df$comparison)
p_cor_df$diff <- factor(p_cor_df$diff, levels = c("up", "dn"))

p_cor_df$label <- ifelse(p_cor_df$corr > 0.5, p_cor_df$RBP, NA)

pdf(file = paste0(plot_out_dir, "boxplot_pvalue_spearman_correlation.pdf"), 6,4)
ggplot(p_cor_df, aes(x = comparison, y = corr, fill = diff))+
  geom_boxplot(alpha = 0.8, position = "dodge2")+
  # geom_violin(alpha = 0.5)+ 
  # geom_dotplot(binaxis = "y",
  #              stackdir = "center",
  #              dotsize = 0.5) +
  
  geom_point(pch = 21, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2)) +
  ylab("Spearman's correlation coefficient") +
  # xlab("Targets of expression quantification") +
  # scale_x_discrete(labels = c("day1_background_day0" = expression("day 0" %->% "day 1"), 
  #                             "day2_background_day1" = expression("day 1" %->% "day 2"), 
  #                
  #                             "day3_background_day2" = expression("day 2" %->% "day 3"), 
  #                             "day5_background_day3" = expression("day 3" %->% "day 5"))) + 
  scale_x_discrete(labels = c("day1_background_day0" = "day 1 vs day 0", 
                              "day2_background_day1" = "day 2 vs day 1", 
                              
                              "day3_background_day2" = "day 3 vs day 2", 
                              "day5_background_day3" = "day 5 vs day 3")) + 
  scale_fill_discrete(name = "Differential exon skipping",
                      labels = c("dn" = "Skipped", 
                              "up" = "Included")) + 
  theme_classic() +
  theme(legend.position = "top", text = element_text(size = 15),
        axis.title.x = element_blank())
dev.off()




#heatmap ####
data_new <- data[df$RBP[which(cor_list >0.5)], ]
rownames(data_new) <- df$RBP[which(cor_list >0.5)]

annot_col <- data.frame(do.call(rbind, strsplit(colnames(data), "_")))
annot_col_new <- annot_col[, -c(1:4), drop = F]
colnames(annot_col_new) <- "RBP binding site"
rownames(annot_col_new) <- colnames(data)
annot_col_new[annot_col_new[,1] == "motif", 1] <- "Consensus motif"
annot_col_new[annot_col_new[,1] == "clip", 1] <- "CLIP-seq"

breaksList = seq(0, 6, by = 0.01)
colors <- colorRampPalette( brewer.pal(9, "Reds") )(length(breaksList))

cols <- c("#F8766D", "#619CFF")
names(cols) <- factor(c("Consensus motif","CLIP-seq"))

anno_colors <- list("RBP binding site" = cols)

pheatmap(data_new,
         show_rownames=T,cluster_rows=F, 
         show_colnames=F,cluster_cols=F,
         annotation_colors=anno_colors,
         # annotation_row=kclasses, annotation_names_row = F,
         annotation_col=annot_col_new, annotation_names_col = F,
         # labels_col = c(expression("day 0" %->% "day 1"), expression("day 1" %->% "day 2"), 
         #                expression("day 2" %->% "day 3"), expression("day 3" %->% "day 5")), 
         angle_col = '45',
         color = colors,
         # annotation_legend = T, 
         fontsize = 14, border_color = NA, 
         breaks = breaksList,
         # gaps_row = head(as.numeric(cumsum(rev(table(kclasses$group)))), -1),
         gaps_col = head(as.numeric(cumsum(rev(table(annot_col$X4)))), -1),
         # filename = paste0(plot_dir_out, 'heatmap_rmaps_pval.without_col_label.pdf'), width = 8, height = 12
)

#balloon plot ####
motif_data <- as.matrix(motif_df %>% select(contains("smallest")))
rownames(motif_data) <- motif_df$RBP

breaksList = seq(0, 6, by = 0.01)
if (d == "up") {
  cols <- colorRampPalette( brewer.pal(9, "Reds") )(length(breaksList))
  plot_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Plots/SE_included_", comp, ".motif_clip.balloon.pdf")
}else{
  cols <- colorRampPalette( brewer.pal(9, "Blues") )(length(breaksList))
  plot_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Plots/SE_skipped_", comp, ".motif_clip.balloon.pdf")
  
}


pdf(plot_fn, 7,12)
ggballoonplot(motif_data, fill = "value",
              size = 8) +
  # gradient_fill(c("blue", "white", "red")) + 
  scale_fill_gradientn(name = "Consensus motif enrichment",
                       colours=cols,
                       # breaks=c(0,6),
                       limits=c(0,6),
                       na.value=cols[length(cols)]
                       ) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.ticks = element_blank(),
        text = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13)) +
  geom_tile(color = "transparent", fill = "transparent") +
  geom_point(data = clip_signif, mapping = aes(x = pos, y = rbp, color = "p < 0.05"),
             shape = 23, fill = "white", size = 1.5, stroke = 1.5) +
  scale_color_manual(name = "CLIP-seq binding site enrichment",
                       values = c("p < 0.05" = "black")) +
  # labs(fill = , color = ) +
  coord_equal()
dev.off()


test <- reshape2::melt(motif_data)

if (d == "up") {
  cols <- colorRampPalette( brewer.pal(9, "Reds") )(length(breaksList))
  plot_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Plots/SE_included_", comp, ".motif_clip.pdf")
}else{
  cols <- colorRampPalette( brewer.pal(9, "Blues") )(length(breaksList))
  plot_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Plots/SE_skipped_", comp, ".motif_clip.pdf")
  
}
pdf(plot_fn, 7,6.8)
ggplot(data=test, aes(x=Var2, y=Var1)) +
  geom_tile(aes(fill = value), width=0.9, height=0.9) +
  geom_point(data = clip_signif, mapping = aes(x = pos, y = rbp, color = "p < 0.05"),
             shape = 23, fill = "white", size = 1.5, stroke = 1.5) +
  scale_fill_gradientn(name = "Consensus motif enrichment",
                       colours=cols,
                       # breaks=c(0,6),
                       limits=c(0,6),
                       na.value = cols[length(cols)]
  ) +
  scale_color_manual(name = "CLIP-seq binding site enrichment",
                     values = c("p < 0.05" = "black")) +
  scale_y_discrete(limits=rev)+ 
  theme_classic() + 
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        plot.title.position = "plot",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.ticks = element_blank(),
        text = element_text(colour = "black", size = 13),
        axis.text.y = element_text(colour = "black", size = 13))
dev.off()







#create a table indicating if clip p < 0.05 at the region with smallest motif p ####
threshold_log <- -log10(0.05)
threshold <- 0.001

comparison <- c("day1_background_day0", "day2_background_day1", "day3_background_day2", "day5_background_day3")
diff <- c("up", "dn")

comp_list <- c()
d_list <- c()
pos_list <- c()
rbp_list <- c()

for (comp in comparison) {
  for (d in diff) {
    motif_fn <- paste0("/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitches_RBPs/rMAPS/Results/", comp, "/SE/pVal.", d, ".vs.bg.RNAmap.txt")
    motif_df <- read.table(motif_fn, sep = '\t', header = T)
    
    motif_df <- motif_df %>%
      filter(rowSums(select(., contains('smallest_p_in_')) < threshold) > 0)
    motif_df <- motif_df %>%
      mutate(RBP = sub("\\..*", "", RBP)) %>%
      arrange(RBP) %>%
      group_by(RBP) %>%
      summarise_all(.funs = list(min))
    
    signif_rbp <- motif_df$RBP
    
    clip_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_pvalue_table/", comp, "/RBP.pVal.", d, ".vs.bg.txt")
    clip_df <- read.table(clip_fn, sep = '\t', header = T)
    
    meta_fn <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Data/RBP_peak_meta.txt"
    meta_df <- read.table(meta_fn, sep = '\t', header = F)
    
    clip_df <- clip_df %>% left_join(meta_df, by = c("RBP" = "V2"))
    
    rownames(clip_df) <- clip_df$V1
    clip_data <- as.matrix(clip_df %>% select(contains("smallest")))
    
    signif_rbp <- intersect(signif_rbp, rownames(clip_df))
    
    
    motif_data <- as.matrix(motif_df %>% filter(RBP %in% signif_rbp) %>% select(contains("smallest")))
    rownames(motif_data) <- signif_rbp
    
    if (nrow(motif_data) > 0) {
      
      for (i in 1:nrow(motif_data)) {
        rbp_oi <- rownames(motif_data)[i]
        region_oi <- colnames(motif_data)[which.min(motif_data[i, ])]
        
        if (clip_data[rbp_oi, region_oi] > threshold_log) {
          comp_list <- c(comp_list, comp)
          d_list <- c(d_list, d)
          pos_list <- c(pos_list, region_oi)
          rbp_list <- c(rbp_list, rbp_oi)
        }
      }
    }
    
    
    
    
    
  }
}

tmp_df <- data.frame("comparison" = comp_list,
                     "diff" = d_list,
                     "pos" = pos_list,
                     "rbp" = rbp_list)
write.table(tmp_df, "/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/cooccurrence_pvalue_motif_clip.txt",
            row.names = F, quote= F, sep = '\t')
