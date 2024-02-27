library(dplyr)
library(ggplot2)
# library(ggpubr)
# library(ggh4x)
library(optparse)

option_list = list(
  make_option(c("-p", "--protein"), type="character", default=NULL, 
              help="rbp name", metavar="character"),
  make_option(c("-c", "--comparison"), type="character", default=NULL, 
              help="day-to-day comparison", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


if (is.null(opt$protein)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (RBP).", call.=FALSE)
}else if (is.null(opt$comparison)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (day-to-day comparison).", call.=FALSE)
}else{
  cat("RBP: ", opt$protein, '\nComparison: ', opt$comparison)
}
  


###rbp-maps outfile metadata

#sum(col(norm_matrix))/ncol = hist (i.o.w., peak density)

#sum(row(norm_matrix)) = coverage

###create peak density and p-value plot

rbp <- opt$protein
# rbp <- "SRSF1"
comp <- opt$comparison
# comp <- "day1_background_day0"
ind <- paste0(rbp, "_", comp)

seq1 <- seq(-50, 250)
seq1 <- seq1[seq1 !=0]

seq2 <- seq(-250, 50)
seq2 <- seq2[seq2 !=0]

res_dir <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Results/CLIPdb_RBP/", rbp, '/', comp, '/')

##up data ####
up_pvalue_fn <- paste0(res_dir, ind, ".up.SE.MATS.JC.txt.pvalues.txt")
up_pvalue_df <- read.table(up_pvalue_fn, sep = '\t', header = F)
colnames(up_pvalue_df) <- c("pos", "pvalue")
up_pvalue_df$pos <- c(seq1, seq2, seq1, seq2)

up_density_fn <- paste0(res_dir, ind, ".up.SE.MATS.JC.txt.hist.txt")
up_density_df <- read.table(up_density_fn, sep = ',', header = F)
colnames(up_density_df) <- "density"



up_df <- cbind(up_pvalue_df, up_density_df)


up_df <- up_df %>% mutate(lineplot_block = c(rep(1, 300), rep(2, 300), rep(3, 300), rep(4,300)))
up_df <- up_df %>% mutate(comparison = "up")


##dn data ####
dn_pvalue_fn <- paste0(res_dir, ind, ".dn.SE.MATS.JC.txt.pvalues.txt")
dn_pvalue_df <- read.table(dn_pvalue_fn, sep = '\t', header = F)
colnames(dn_pvalue_df) <- c("pos", "pvalue")
dn_pvalue_df$pos <- c(seq1, seq2, seq1, seq2)

dn_density_fn <- paste0(res_dir, ind, ".dn.SE.MATS.JC.txt.hist.txt")
dn_density_df <- read.table(dn_density_fn, sep = ',', header = F)
colnames(dn_density_df) <- "density"

dn_df <- cbind(dn_pvalue_df, dn_density_df)


dn_df <- dn_df %>% mutate(lineplot_block = c(rep(1, 300), rep(2, 300), rep(3, 300), rep(4,300)))
dn_df <- dn_df %>% mutate(comparison = "dn")

##bg data ####
bg_density_fn <- paste0(res_dir, ind, ".bg.SE.MATS.JC.txt.hist.txt")
bg_density_df <- read.table(bg_density_fn, sep = ',', header = F)
colnames(bg_density_df) <- "density"

bg_df <- bg_density_df
bg_df$pos <- up_df$pos
bg_df$pvalue <- NA

bg_df <- bg_df %>% mutate(lineplot_block = c(rep(1, 300), rep(2, 300), rep(3, 300), rep(4,300)))
bg_df <- bg_df %>% mutate(comparison = "bg")

bg_df <- bg_df[, colnames(up_df)]


##combined into one df
df <- rbind(up_df, dn_df, bg_df)
cols <- c("up" = "red", "dn" = "blue", "bg" = "black")


coeff <- ceiling((max(df$pvalue, na.rm = T)/max(df$density, na.rm = T))/10)*10
plot_fn <- paste0("/project/Neurodifferentiation_System/owlmayerTemporary/derek/CLiP_rMAPs/Plots/ClIPdb_RBP/", ind, ".pdf")
pdf(plot_fn, 11.2,2.9)
ggplot(df, aes(x = pos, color = comparison)) +
  geom_line( aes(y=density), linewidth = 0.7, linetype = "solid") + 
  geom_line( aes(y=pvalue/coeff), linewidth = 0.3, linetype = "dashed") + 
  geom_vline(xintercept = 0, linewidth = 0.1) +
  scale_colour_manual(values = cols) + 
  scale_y_continuous(
    # limits = c(0, 0.11),
    # Features of the first axis
    name = "Peak density",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="-log10(p-value)")
  ) + 
  facet_wrap(~lineplot_block, nrow = 1, scales = "free_x") +
  scale_x_continuous(
    breaks = c(-250, -125, -50, 0, 50, 125, 250),
    expand = c(0, 0)
  )+
  theme_bw() +
  theme( strip.background = element_blank(),
         strip.text = element_blank(),
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
         panel.spacing.x = unit(8, "mm"),
         plot.title = element_text(hjust = 0.5),
         text = element_text(size = 15, colour = "black"),
         axis.text = element_text(size = 12, colour = "black"),
         axis.ticks = element_blank(),
         axis.title.x=element_blank()
         
  )+
  ggtitle(paste0(rbp, ": ", comp)) +
  guides(colour=FALSE) 


dev.off()

