#loading the packages
library(tidyverse)
library(tximport)
library(DESeq2)
library(ashr)
library(annotables)
library(ape)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(scico)
library(colorspace)
library(RColorBrewer)
library(pheatmap)
library(conflicted)
library(forcats)
library(ggfortify)
library(ggpubr)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("slice", "dplyr")

### Prepare data ####
#reading in the paths and metadata ####
annotfn <- '/project/PausingDynamics/GeneralResources/Genecode29/gencode.v29.annotation.sorted.gff3'
exp_dir<-'/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/'
brain_dir <- '/project/Neurodifferentiation_System/Published_data/BrainAtlas/'
res_dir_out<-'/project/Neurodifferentiation_System/owlmayerTemporary/derek/atlas/Results/'
plot_dir_out<-'/project/Neurodifferentiation_System/owlmayerTemporary/derek/atlas/Plots/'

samples <- c('OJ63','OJ64','OJ65','OJ66','OJ67','OJ68','OJ69','OJ70','OJ71','OJ72','OJ73','OJ74',
             'OJ75','OJ76','OJ77','OJ79','OJ80','OJ81','OJ83','OJ85','OJ86','OJ87','OJ88','OJ89')
cond <- c('day0','day0','day0','day0','day0','day1','day1','day1','day1','day1','day2','day2',
          'day2','day2','day2','day3','day3','day3','day3','day5','day5','day5','day5','day5')
rep <- c('rep1','rep2','rep3','rep4','rep5','rep1','rep2','rep3','rep4','rep5','rep1','rep2',
          'rep3','rep4','rep5','rep1','rep2','rep3','rep5','rep1','rep2','rep4','rep5','rep6')
fns <- paste0(exp_dir,samples,'.genes.results')
pval <- 0.05
l2fc <- log2(1.5)

qvalgo<-0.05

breaksList = seq(0.45, 0.7, by = 0.01)
breaksList = seq(0.4, 0.8, by = 0.01)
breaksList = seq(0.4, 0.7, by = 0.02)
breaksList = seq(-2, 3, by = 0.1)
colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
colors <- colorRampPalette( brewer.pal(9, "Reds") )(length(breaksList))


#loading the Brain Atlas data ####
brain_exp <- read.csv(paste0(brain_dir, 'expression_matrix.csv'), header = F, row.names = 1)
brain_sample <- read.csv(paste0(brain_dir, 'columns_metadata.csv'))
brain_gene <- read.csv(paste0(brain_dir, 'rows_metadata.csv'))

rownames(brain_exp) <- brain_gene$ensembl_gene_id

#filter out unnecessary structures and ages ####
#structure metadata: retain freq>=24 structures -> top 16
# ggplot(brain_sample, aes(x = fct_infreq(structure_acronym), fill = age))+
  # geom_bar(stat = "count", position = "stack")

struc_meta <- matrix(nrow = length(unique(brain_sample$structure_acronym)), ncol = 3)
struc_meta <- as.data.frame(struc_meta)
colnames(struc_meta) <- c("structure", "index", "freq")

struc_meta$structure <- unique(brain_sample$structure_acronym)
struc_meta$index <- lapply(struc_meta$structure, function(x){y = which(brain_sample$structure_acronym == x)})
struc_meta$freq <- unlist(lapply(struc_meta$index, function(x){y = length(x)}))


# pdf(paste0(plot_dir_out, "brain_struc_freq.pdf"), 12, 4)
# ggplot(struc_meta, aes(x = reorder(structure, unlist(freq), decreasing = T), y = freq)) +
#   geom_bar(stat = "identity") +
#   xlab("brain structure") +
#   ylab("frequency")
# dev.off()

struc_meta <- struc_meta[struc_meta$freq >= 24, ]
struc_meta$full_name <- lapply(struc_meta$structure, function(x){y = brain_sample[brain_sample$structure_acronym == x, ]$structure_name[1]})


#age metadata: retain num_structure>1 then number of ages -> 30
age_meta <- matrix(nrow = length(unique(brain_sample$age)), ncol = 2)
age_meta <- as.data.frame(age_meta)
colnames(age_meta) <- c("age", "index")
age_meta$age <- unique(brain_sample$age)
age_meta$index <- lapply(age_meta$age, function(x){y = which(brain_sample$age == x)})
age_meta$num_struc <- unlist(lapply(age_meta$index, function(x){y = length(unique(brain_sample$structure_acronym[unlist(x)]))}))

age_meta <- age_meta[age_meta$num_struc > 1, ]


# brain_of_interest <- intersect(unlist(age_meta$index), unlist(struc_meta$index))


#group: age + structure
# brain_sample$age <- sapply(brain_sample$age, function(x){y=gsub(" ", "", x)})
# brain_sample$group <- paste(brain_sample$age, brain_sample$structure_acronym, sep = "_")
# 
# colnames(brain_exp) <- brain_sample$group
# 
# unique(brain_sample$age)
# unique(brain_sample$structure_acronym)
# unique(brain_sample$group)
# 
# brain_sample$group[duplicated(brain_sample$group)]



#convert RPKM to TPM ####
brain_exp <- apply(brain_exp, MARGIN = 2, function(x){y = x/sum(x)*10^6})
brain_exp <- as.data.frame(brain_exp)


# 
# age_subgroup <- which(brain_sample$age %in% unique(brain_sample$age)[1] == T)
# corMatrix <- cor(brain_exp[, subgroup])
#w.r.t structure
# rownames(corMatrix) <- brain_sample$structure_acronym[subgroup]
# colnames(corMatrix) <- brain_sample$structure_acronym[subgroup]
# colors <- colorRampPalette( brewer.pal(9, "Reds") )(255)
# pheatmap(corMatrix,col=colors)




#loading the RSEM results ####
txi<-tximport(fns, type="rsem", txIn = FALSE, txOut = FALSE)
txi$length[txi$length == 0] <- 1


#load TPM quantification of cell expression ####
cell_exp_mat <- txi$abundance
colnames(cell_exp_mat)<-paste(cond,rep,sep='_')
cell_exp_mat<- as.data.frame(cell_exp_mat)

cell_exp_mat <- cell_exp_mat[-which(grepl("PAR_Y", rownames(cell_exp_mat), fixed = T)), ]


#uniform format of ensembl gene id with Brain Atlas data ####
for (i in 1:nrow(cell_exp_mat)) {
  rownames(cell_exp_mat)[i] <- unlist(strsplit(rownames(cell_exp_mat)[i], '\\.'))[1]
}


#common set of genes - 45603
common_gene <- base::intersect(rownames(cell_exp_mat), rownames(brain_exp))
brain_exp <- brain_exp[common_gene, ]
cell_exp_mat <- cell_exp_mat[common_gene, ]

#QC of our cell samples and brain atlas samples
# test <- t(as.matrix(cell_exp_mat))
# test <- as.data.frame(test)
# pca_res <- prcomp(test)
# test <- test %>% mutate(day = unlist(lapply(rownames(test), function(x) {y<- substr(x, 1, 4)})))
# 
# autoplot(pca_res, data = test, colour = "day")
# 
# tmp <- t(as.matrix(brain_exp))
# tmp <- as.data.frame(tmp)
# 
# pca_res <- prcomp(tmp)
# 
# tmp <- tmp %>% mutate(age = factor(brain_sample$age, levels = age_meta$age))
# tmp <- tmp %>% mutate(struc = brain_sample$structure_acronym)
# pdf(paste0(plot_dir_out, 'brain_pca_age.pdf'), 8, 6)
# autoplot(pca_res, data = tmp, colour = "age")
# dev.off()
# pdf(paste0(plot_dir_out, 'brain_pca_structure.pdf'), 8, 6)
# autoplot(pca_res, data = tmp, colour = "struc")
# dev.off()




#distribution of TPM quantification - extremely right-skewed
# brain_exp_mat <- as.matrix(brain_exp)
# plot(density(log(c(brain_exp_mat)+1)))
# max(brain_exp)
# 
# cell_exp_mat <- as.matrix(cell_exp)
# plot(density(log(c(cell_exp_mat)+1)))
# max(cell_exp)


#generate filtered gene list ####
#filter out low-expression genes with mean/median TPM <= 1 
filt_common_gene <- c()
for (i in 1:length(common_gene)) {
  v <- cell_exp_mat[i,]
  u <- brain_exp[i,]
  names(v) <- NULL
  v <- unlist(c(v))
  names(u) <- NULL
  u <- unlist(c(u))
  # if (mean(v, na.rm = T) > 1 & mean(u, na.rm = T) > 1) {
  if (median(v, na.rm = T) > 1 & median(u, na.rm = T) > 1) {
  #if (mean(v) > 1) {
    filt_common_gene <- c(filt_common_gene, common_gene[i])
  }
}

# cell_exp_mat <- as.matrix(cell_exp[filt_common_gene, ])
# plot(density(log(c(cell_exp_mat)+1)))



### Define functions ####

#define a function for sample-wise correlation calculation ####
calCorMat <- function(row_expression_mat = cell_exp_mat, col_expression_mat = brain_exp,
                      cor_method = "pearson", log = F, gene_list = c(1:length(common_gene)),
                      day = "day5") {
  row_expression_mat <- row_expression_mat %>% select(contains(day))
  
  cor_mat <- matrix(nrow = ncol(row_expression_mat), ncol = ncol(col_expression_mat))
  
  if (!log) {
    for (i in (1:ncol(cor_mat))) {
      for (j in (1:nrow(cor_mat))) {
          cor_mat[j,i] <- cor(row_expression_mat[gene_list, j], col_expression_mat[gene_list,i], method = cor_method)
      }
    }
  }else{
    for (i in (1:ncol(cor_mat))) {
      for (j in (1:nrow(cor_mat))) {
        cor_mat[j,i] <- cor(log2(row_expression_mat[gene_list, j]+1), log2(col_expression_mat[gene_list,i]+1), method = cor_method)
      }
    }
  }
  print(paste0('Max corr: ', max(cor_mat, na.rm = T)))
  print(paste0('Min corr: ', min(cor_mat, na.rm = T)))
  
  cor_mat
}


#define a function for heatmap matrix based on sample-wise correlation ####
calHeatmapMat <- function(row_meta = age_meta, col_meta = struc_meta, 
                          corr_mat = cor_mat, replace_NA = F, z_norm = T){
  
  heatmap_mat <- matrix(nrow = nrow(row_meta), ncol = nrow(col_meta))
  
  for (i in 1:nrow(heatmap_mat)) {
    for (j in 1:ncol(heatmap_mat)) {
      #take the median to turn sample-to-sample correlation into group-to-group correlation
      group_index <- base::intersect(unlist(row_meta$index[i]), unlist(col_meta$index[j]))
      heatmap_mat[i,j] <- median(corr_mat[, group_index])
    }
  }
  
  if(replace_NA) {
    for (i in 1:nrow(heatmap_mat)) {
      for (j in 1:ncol(heatmap_mat)) {
        if (is.nan(heatmap_mat[i,j])) {
          heatmap_mat[i,j] <- median(heatmap_mat[i,], na.rm = T)
        }
      }
    }
  }
  
  if(z_norm) {
    #row-wise norm
    # for (i in 1:nrow(heatmap_mat)) {
    #   heatmap_mat[i, ] <- (heatmap_mat[i,]-mean(heatmap_mat[i,], na.rm = T))/sd(heatmap_mat[i,], na.rm = T)
    # }
    
    #matrix-wise norm
    print(paste0('Max corr: ', max(heatmap_mat, na.rm = T)))
    print(paste0('Min corr: ', min(heatmap_mat, na.rm = T)))
    print(paste0('mean: ', mean(heatmap_mat, na.rm = T)))
    print(paste0('sd: ', sd(heatmap_mat, na.rm = T)))
    
    tmp <- heatmap_mat
    for (i in 1:nrow(heatmap_mat)) {
      for (j in 1:ncol(heatmap_mat)) {
        heatmap_mat[i,j] <- (tmp[i,j]-mean(tmp, na.rm = T))/sd(tmp, na.rm = T)
      
      }
    }
    print(paste0('Max scaled value: ', max(heatmap_mat, na.rm = T)))
    print(paste0('Min scaled value: ', min(heatmap_mat, na.rm = T)))
    
    
  }else{
    print(paste0('Max corr: ', max(heatmap_mat, na.rm = T)))
    print(paste0('Min corr: ', min(heatmap_mat, na.rm = T)))
    
  }
  
  heatmap_mat
  
}

#define a function for comparison of correlation coefficients between two days  ####
compareCor <- function(cor_mat_1, day1 = 'day0', cor_mat_2, day2 = 'day5', age_info = age_meta, struc_info = struc_meta, sample_info = brain_sample) {
  
  valid_comb <- c()
  for (i in 1:nrow(age_info)) {
    for (j in 1:nrow(struc_info)) {
      valid_comb <- c(valid_comb, intersect(unlist(age_info$index[i]), unlist(struc_info$index[j])))
    }
  }
  valid_comb <- sort(valid_comb)
  
  day_list <- c(rep(day1, length(valid_comb)*nrow(cor_mat_1)), rep(day2, length(valid_comb)*nrow(cor_mat_2)))
  
  age_list <- c()
  for (i in 1:length(valid_comb)) {
    age <- sample_info$age[valid_comb[i]]
    tmp <- rep(age, nrow(cor_mat_1))
    age_list <- c(age_list, tmp)
  }
  for (i in 1:length(valid_comb)) {
    age <- sample_info$age[valid_comb[i]]
    tmp <- rep(age, nrow(cor_mat_2))
    age_list <- c(age_list, tmp)
  }
  
  cor_list <- c(as.vector(cor_mat_1[, valid_comb]), as.vector(cor_mat_2[, valid_comb]))

  
  df <- data.frame(day = day_list, age = age_list, cor = cor_list)
  
  df
  
}

###Compare correlation coefficients between day0 and day5 ####
cor_mat_0 <- calCorMat(cor_method = "spearman", gene_list = filt_common_gene, day = 'day0')
cor_mat_5 <- calCorMat(cor_method = "spearman", gene_list = filt_common_gene, day = 'day5')

df <- compareCor(cor_mat_1 = cor_mat_0, cor_mat_2 = cor_mat_5)

pdf(paste0(plot_dir_out, 'corr_compare_0_vs_5.pdf'), 15, 5)
ggplot(df, aes(x = factor(age, age_meta$age), y = cor, fill = day), show.legend = FALSE) +
  geom_boxplot() +
  stat_compare_means(aes(group = day), method = "wilcox.test", label =  "p.signif") +
  theme_bw() +
  xlab("Age of Brain Atlas samples") +
  ylab("Spearman's correlation") +
  guides(fill=guide_legend(title="Day"))
dev.off()


###Plot heatmap for individual day ####
day <- "day5"
##all 45603 genes ####

#sample-wise Pearson's cor on raw TPM - max. 0.14, min. 0.004 ####
cor_mat <- calCorMat()

heatmap_mat <- calHeatmapMat()

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = T, 
         col=colors, main = "sample-wise Pearson on raw TPM of all genes")



#sample-wise Spearman's cor on Raw TPM - max. 0.72, min. 0.58 ####
cor_mat <- calCorMat(cor_method = "spearman", day = day)

#calculate within-group variance to support the selection of median
# group <- c()
# for ( i in 1:nrow(age_meta) ){
#   for (j in 1:nrow(struc_meta)) {
#     group <- c(group, paste(age_meta$age[i], struc_meta$structure[j], sep = "_"))
#   }
# }
# 
# group_var <- matrix(nrow = length(group), ncol = 2)
# group_var[, 1] <- group
# for (i in 1:nrow(age_meta)) {
#   for (j in 1:nrow(struc_meta)) {
#     #take the median to turn sample-to-sample correlation into group-to-group correlation
#     group_index <- intersect(unlist(age_meta$index[i]), unlist(struc_meta$index[j]))
#     if(i == 25 & j == 1) {print(var(as.vector(cor_mat[, group_index])))}
#     if(length(group_index) == 0) {
#       group_var[(i-1)*j+j,2] <- NA
#     }else{
#       group_var[(i-1)*nrow(struc_meta)+j,2] <- var(as.vector(cor_mat[, group_index]))
#     }
#     
#   }
# }
# 
# group_var <- as.data.frame(group_var)
# colnames(group_var) <- c("group", "variance")
# group_var$variance <- as.numeric(group_var$variance)
# 
# pdf(paste0(plot_dir_out,'within-group_variance.pdf'), 4, 3)
# ggplot(group_var, aes(x = variance)) +
#   geom_density(na.rm = T) +
#   xlab("within-group variance")
# dev.off()


heatmap_mat <- calHeatmapMat(z_norm = F)
median(heatmap_mat, na.rm = T)
max(heatmap_mat, na.rm = T)
min(heatmap_mat, na.rm = T)
rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


# heatmap_mat_nonZ <- calHeatmapMat(z_norm = F)
# rownames(heatmap_mat_nonZ) <- age_meta$age
# colnames(heatmap_mat_nonZ) <- struc_meta$structure
# heatmap_mat_nonZ <- as.data.frame(heatmap_mat_nonZ)
# write.csv(heatmap_mat_nonZ, paste0(res_dir_out, "day5_spearman_all_raw.csv"), row.names = T)

red_breaksList = seq(median(heatmap_mat, na.rm = T), max(heatmap_mat, na.rm = T), by = 0.002)
blue_breaksList = seq(min(heatmap_mat, na.rm = T), median(heatmap_mat, na.rm = T), by = 0.002)


#median as the separation line
colors <- c(colorRampPalette( rev(brewer.pal(9, "Blues")) )(length(blue_breaksList)), 
            colorRampPalette( brewer.pal(9, "Reds") )(length(red_breaksList)))



pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors,
         fontsize = 8, 
         # breaks = breaksList,
         filename=paste0(plot_dir_out, 'heatmap_spearman_all_raw_', day, '.pdf'), width=4, height=6
         )


#sample-wise Pearson's cor on log2 TPM - max. 0.72, min. 0.54 ####
cor_mat <- calCorMat(log = T, day = 'day5')

heatmap_mat <- calHeatmapMat(z_norm = F)

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors, fontsize = 8, breaks = breaksList,
         filename=paste0(plot_dir_out, 'heatmap_pearson_all_log_raw_', day, '.pdf'), width=4, height=6)


##filtered 14345 genes ####
#stats (mean & sd) to scale over all days ####
# days <- c('day0', 'day1', 'day2', 'day3', 'day5')
days <- c('day0', 'day5')

heatmap_mat_list <- c()
cor_value_list <- c()

for(i in 1:length(days)) {
  cor_mat <- calCorMat(cor_method = "spearman", gene_list = filt_common_gene,
                       day = days[i])
  heatmap_mat <- calHeatmapMat(z_norm = F)
  heatmap_mat_list <- c(heatmap_mat_list, list(heatmap_mat))
  cor_value_list <- c(cor_value_list, heatmap_mat)
}

stats <- c(mean(cor_value_list, na.rm = T), sd(cor_value_list, na.rm = T))

for (i in 1:length(days)) {
  heatmap_mat <- heatmap_mat_list[[i]]
  for (j in 1:nrow(heatmap_mat)) {
    for (k in 1:ncol(heatmap_mat)) {
      heatmap_mat[j,k] <- (heatmap_mat[j,k]-stats[1])/stats[2]
      
    }
  }
  print(paste0('Max scaled value: ', max(heatmap_mat, na.rm = T)))
  print(paste0('Min scaled value: ', min(heatmap_mat, na.rm = T)))
  rownames(heatmap_mat) <- age_meta$age
  colnames(heatmap_mat) <- struc_meta$structure
  heatmap_mat_list[[i]] <- heatmap_mat
}


pheatmap(heatmap_mat_list[[1]], na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors,
         fontsize = 8, breaks = breaksList,
         filename=paste0(plot_dir_out, paste0('heatmap_spearman_filtered_overall.scaled_day0.pdf')), width=4, height=6
        
)
pheatmap(heatmap_mat_list[[2]], na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors,
         fontsize = 8, breaks = breaksList,
         filename=paste0(plot_dir_out, paste0('heatmap_spearman_filtered_overall.scaled_day5.pdf')), width=4, height=6
         
)

#
#sample-wise Pearson's cor on raw TPM - max. , min.  ####
cor_mat <- calCorMat(gene_list = filt_common_gene, day = 'day5')

heatmap_mat <- calHeatmapMat(z_norm = F)

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors, fontsize = 8,
         filename=paste0(plot_dir_out, paste0('heatmap_pearson_filtered_raw_', day, '.pdf')), width=4, height=6
)         

#sample-wise Spearman's cor on Raw TPM - max. , min.  ####
cor_mat <- calCorMat(cor_method = "spearman", gene_list = filt_common_gene, day = day)

#calculate within-group variance to support the selection of median

# group <- c()
# for ( i in 1:nrow(age_meta) ){
#   for (j in 1:nrow(struc_meta)) {
#     group <- c(group, paste(age_meta$age[i], struc_meta$structure[j], sep = "_"))
#   }
# }
# 
# group_var <- matrix(nrow = length(group), ncol = 2)
# group_var[, 1] <- group
# for (i in 1:nrow(age_meta)) {
#   for (j in 1:nrow(struc_meta)) {
#     #take the median to turn sample-to-sample correlation into group-to-group correlation
#     group_index <- intersect(unlist(age_meta$index[i]), unlist(struc_meta$index[j]))
#     if(i == 25 & j == 1) {print(var(as.vector(cor_mat[, group_index])))}
#     if(length(group_index) == 0) {
#       group_var[(i-1)*j+j,2] <- NA
#     }else{
#       group_var[(i-1)*nrow(struc_meta)+j,2] <- var(as.vector(cor_mat[, group_index]))
#     }
#     
#   }
# }
# 
# group_var <- as.data.frame(group_var)
# colnames(group_var) <- c("group", "variance")
# group_var$variance <- as.numeric(group_var$variance)
# 
# pdf(paste0(plot_dir_out,'within-group_variance.pdf'), 4, 3)
# ggplot(group_var, aes(x = variance)) +
#   geom_density(na.rm = T) +
#   xlab("within-group variance")
# dev.off()


heatmap_mat <- calHeatmapMat(z_norm = F)

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure
# which(heatmap_mat == max(heatmap_mat, na.rm = T), arr.ind = T)


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors,
         fontsize = 8, breaks = breaksList,
         # filename=paste0(plot_dir_out, paste0('heatmap_spearman_filtered_raw_', day, '.pdf')), width=4, height=6
         # filename=paste0(plot_dir_out, paste0('heatmap_spearman_filtered_scaled_', day, '.pdf')), width=4, height=6
         
         )

#
#sample-wise Pearson's cor on log2 TPM - max. , min.  ####
cor_mat <- calCorMat(log = T, gene_list = filt_common_gene, day = 'day5')

heatmap_mat <- calHeatmapMat(z_norm = T)

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors, fontsize = 8, #breaks = breaksList,
         #filename=paste0(plot_dir_out, paste0('heatmap_pearson_filtered_log_raw_', day, '.pdf')), width=4, height=6
)










##top 500/1000/2000 up-regulated genes for day5-vs-day0 ####

DEfn='/project/Neurodifferentiation_System/Analysis_NGN3/DEseq_genes/Results/day5_background_day0/diff_expression.csv'
de <- read.csv(DEfn, header = T, row.names = 1)


de_gene_list <- de$ensgene
# de_gene_list <- intersect(de_gene_list, common_gene)
de_gene_list <- intersect(de_gene_list, filt_common_gene)

up_de <- de[de$log2FoldChange > 0, ]

up_gene_list <- up_de[order(up_de$padj), ]$ensgene

up_gene_list <- intersect(up_gene_list, filt_common_gene)


up500_gene_list <- intersect(up_gene_list, filt_common_gene)[1:500]
up1000_gene_list <- intersect(up_gene_list, filt_common_gene)[1:1000]
up2000_gene_list <- intersect(up_gene_list, filt_common_gene)[1:2000]



#comparison of correlation coefficient distribution ####
compareCCDistribution <- function(GOI_list = up_gene_list, 
                                  number_list = c(500, 1000, 2000), 
                                  cor_method = "pearson", log = F) {
  value = c()
  gene_list = c()
  
  cor_mat <- calCorMat(cor_method = cor_method, log = log)
  cor <- as.vector(cor_mat)
  value = c(value, cor)
  gene_list = c(gene_list, rep("all genes (45603)", length(cor)))
  cor_mat <- calCorMat(cor_method = cor_method, log = log, gene_list = filt_common_gene)
  cor <- as.vector(cor_mat)
  value = c(value, cor)
  gene_list = c(gene_list, rep("filtered genes (17759)", length(cor)))
  for (i in 1:length(number_list)) {
    n = number_list[i]
    genes = GOI_list[1:n]
    cor_mat <- calCorMat(cor_method = cor_method, log = log, gene_list = genes)
    cor <- as.vector(cor_mat)
    value = c(value, cor)
    gene_list = c(gene_list, rep(paste0("top ", n, " upreg 5vs0-upreg-genes"), length(cor)))
  }
  df <- data.frame(value, gene_list)
  
  df
}

df <- compareCCDistribution(cor_method = "spearman")
fn = paste0(plot_dir_out, paste0("corr_dist_diffGenes_spearman.pdf"))
pdf(fn, 7, 5)
ggplot(df, aes(x = value, color = gene_list)) +
  geom_density() +
  xlab("value of correlation coefficient")
dev.off()

df <- compareCCDistribution(cor_method = "pearson", log = T)
fn = paste0(plot_dir_out, paste0("corr_dist_diffGenes_pearson_logTPM.pdf"))
pdf(fn, 7, 5)
ggplot(df, aes(x = value, color = gene_list)) +
  geom_density() +
  xlab("value of correlation coefficient")
dev.off()



#sample-wise Pearson's cor on raw TPM - max. 0.65, min. 0 ####
# cor_mat <- calCorMat(gene_list = up_gene_list)
cor_mat <- calCorMat(gene_list = de_gene_list, day = 'day5')
max(cor_mat)
min(cor_mat)
heatmap_mat <- calHeatmapMat()

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = T, 
         col=colors, main = "sample-wise Pearson on raw TPM of top 500 up-regulated genes")

#sample-wise Spearman's cor on raw TPM - max. 0.36, min. 0.12 ####
# cor_mat <- calCorMat(cor_method = "spearman", gene_list = up_gene_list)
cor_mat <- calCorMat(cor_method = "spearman", gene_list = de_gene_list, day = 'day5')
max(cor_mat)
min(cor_mat)
heatmap_mat <- calHeatmapMat(z_norm = F)

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = F, 
         col=colors, main = "sample-wise Spearman on raw TPM of top 500 up-regulated genes")

#sample-wise Pearson's cor on log TPM - max. 0.37, min. 0.13 ####
# cor_mat <- calCorMat(log = T, gene_list = up_gene_list)
cor_mat <- calCorMat(log = T, gene_list = de_gene_list, day = 'day5')
max(cor_mat)
min(cor_mat)
heatmap_mat <- calHeatmapMat()

rownames(heatmap_mat) <- age_meta$age
colnames(heatmap_mat) <- struc_meta$structure


pheatmap(heatmap_mat, na.color="grey", cluster_row = F, cluster_cols = T, 
         col=colors, main = "sample-wise Pearson on log2 TPM of top 500 up-regulated genes")






