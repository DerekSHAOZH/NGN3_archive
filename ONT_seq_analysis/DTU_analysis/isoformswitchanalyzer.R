library("IsoformSwitchAnalyzeR")
library(ggplot2)
library(dplyr)
library(ggalluvial)
# install.packages("ggthemes")
# library(ggthemes)
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 15
mycolors <- colorRampPalette(brewer.pal(8, "Set3"))(nb.cols)

blue <- colorRampPalette(c('#A6CEE3', "#1F78B4"))(4)
red <- colorRampPalette(c("#B2DF8A", "#33A02C"))(3)
purple <- colorRampPalette(c("#FB9A99", "#E31A1C"))(2)
orange <- colorRampPalette(c("#FDBF6F", "#FF7F00"))(2)
grey <- colorRampPalette(c("#CAB2D6", "#6A3D9A"))(2)

colors <- c(blue, red, purple, orange, grey)

#set seed for reproductive results
set.seed(1) 


cond <- c('day0','day0','day0','day3','day3','day3','day5','day5','day5')
rep <- rep(c(1:3),3)
myDesign <- data.frame(sampleID = paste(cond,rep,sep='_'),condition = cond)

annotf <- "/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/GffCompare/nanopore.combined.filt.gtf"
# annotf <- "/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/GffCompare/nanopore.combined.gtf"
fastaf <- "/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Pinfish/corrected_transcriptome_polished_collapsed_nonredundant_tcons.fas"
tcpf <- "/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_deseq2norm.txt"
direOut <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/"
plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Plots/'

###Data preparation ####
counts<- read.csv(tcpf)
counts<- counts[!duplicated(counts$transcript_id), ]
iso_info <- counts[, 11:ncol(counts)]
tcmf<-counts%>%select(starts_with("day"))%>% select(order(colnames(.)))
tcmf<-as.matrix(tcmf)
rownames(tcmf)<- counts$transcript_id

#number of isoforms: 49557
#number of isoforms (annotated): 26250
length(unique(counts$ref_transcript))
#number of genes (annotated): 13928
length(unique(counts$gene_name))


###Create the object for IsoformSwtichAnalyzeR ####
asswitchlist <- importRdata(isoformExonAnnoation = annotf, isoformNtFasta = fastaf, isoformCountMatrix = tcmf, designMatrix = myDesign)
##########################
#Filter genes: ####
#1) gene TPM > 1 ####
#2) isoform TPM > 0 ####
#3) multi-exon genes  ####
##########################
asswitchlist <- preFilter(switchAnalyzeRlist = asswitchlist,geneExpressionCutoff = 1,isoformExpressionCutoff = 0,removeSingleIsoformGenes = TRUE)


test <- asswitchlist$isoformFeatures
test <- test %>% left_join(iso_info, by = c("isoform_id" = "transcript_id"))
#number of isoforms after filtering: 42028
length(unique(test$isoform_id))
#number of isoforms after filtering (annotated): 21447
length(unique(test$ref_transcript))
#number of multi-exon genes (annotated): 10071
length(unique(test$gene_name.y))


###Differential transcript usage (DTU) analysis ####
# asswitchlist_analyzed <- isoformSwitchTestDRIMSeq(switchAnalyzeRlist = asswitchlist, reduceToSwitchingGenes=FALSE)
asswitchlist_analyzed <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = asswitchlist, reduceToSwitchingGenes=FALSE)

###Analyze ORF ####
#RNA sequence
asswitchlist_analyzed <- analyzeORF(asswitchlist_analyzed, orfMethod = "longest") #, genomeObject = Hsapiens)
# write.csv(asswitchlist_analyzed$orfAnalysis, paste0(direOut,'ORF_analysis.csv'))
#AA sequence
#export the RNA & AA sequences
# asswitchlist_analyzed <- extractSequence(asswitchlist_analyzed, pathToOutput = direOut, writeToFile = TRUE, 
#                                          alsoSplitFastaFile = TRUE, removeLongAAseq= TRUE)
asswitchlist_analyzed <- extractSequence(asswitchlist_analyzed, writeToFile = FALSE, 
                       removeLongAAseq= TRUE)

###Import results from external tools for consequence analysis ####
##CPAT ####
#Turn result files obtained from local execution to the format from web-server
# pathToCPATresultFile <- "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/CPAT/ONT_CP.ORF_prob.best.tsv"
# res <- read.csv(pathToCPATresultFile, sep = "\t")
# colnames(res)[which(colnames(res) == "seq_ID")] <- "Sequence Name"
# colnames(res)[which(colnames(res) == "mRNA")] <- "RNA size"
# colnames(res)[which(colnames(res) == "ORF")] <- "ORF size"
# colnames(res)[which(colnames(res) == "Fickett")] <- "Ficket Score"
# colnames(res)[which(colnames(res) == "Hexamer")] <- "Hexamer Score"
# colnames(res)[which(colnames(res) == "Coding_prob")] <- "Coding Probability"
# res$Index <- c(0:(nrow(res)-1))
# res <- res[, c(12, 1, 3, 8, 9, 10, 11, 2)]
# write.table(res, "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/CPAT/ONT_CP_result.txt", row.names = FALSE, quote = FALSE, sep = "\t")


asswitchlist_analyzed <- analyzeCPAT(
  switchAnalyzeRlist = asswitchlist_analyzed,
  pathToCPATresultFile = "/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/CPAT/ONT_CP_result.txt",
  codingCutoff         = 0.725, # the coding potential cutoff we suggested for human
  removeNoncodinORFs   = TRUE   # because ORF was predicted de novo
)


##Pfam ####
asswitchlist_analyzed <- analyzePFAM(
  switchAnalyzeRlist   = asswitchlist_analyzed,
  pathToPFAMresultFile = '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/Pfam/ONT_domain.txt',
  showProgress=TRUE
)

##SignalP ###
asswitchlist_analyzed <- analyzeSignalP(
  switchAnalyzeRlist       = asswitchlist_analyzed,
  pathToSignalPresultFile  = '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SignalP/ONT_SignalP_summary.signalp5'
)


asswitchlist_analyzed


###Analyze switch consequence ####
consequencesOfInterest <- c("tss", "tts", "last_exon", "isoform_seq_similarity", 'exon_number', 'intron_structure', 'coding_potential',
                            "ORF_seq_similarity", "ORF_length", "5_utr_length", '3_utr_seq_similarity', 'NMD_status',
                            'domains_identified', 'domain_length', 'signal_peptide_identified')

asswitchlist_analyzed_final <- analyzeSwitchConsequences(asswitchlist_analyzed, consequencesToAnalyze=consequencesOfInterest,
                                 removeNonConseqSwitches = FALSE)
#add consequence to intron_structure
asswitchlist_analyzed_final$switchConsequence$switchConsequence[asswitchlist_analyzed_final$switchConsequence$featureCompared == "intron_structure" & asswitchlist_analyzed_final$switchConsequence$isoformsDifferent == TRUE] <- "different intron structure"


conseq_summary <- asswitchlist_analyzed_final$switchConsequence
test <- conseq_summary %>% left_join(iso_info, by = c("isoformUpregulated" = "transcript_id"))
test <- test %>% left_join(iso_info, by = c("isoformDownregulated" = "transcript_id"))
colnames(test)[which(colnames(test) == "gene_name.y")] <- "UP_gene_name"
colnames(test)[which(colnames(test) == "ref_transcript.x")] <- "Up_ref_transcript"
colnames(test)[which(colnames(test) == "gene_id.y")] <- "Up_gene_id"
colnames(test)[which(colnames(test) == "gene_type.x")] <- "Up_gene_type"
colnames(test)[which(colnames(test) == "gene_name")] <- "Down_gene_name"
colnames(test)[which(colnames(test) == "ref_transcript.y")] <- "Down_ref_transcript"
colnames(test)[which(colnames(test) == "gene_id")] <- "Down_gene_id"
colnames(test)[which(colnames(test) == "gene_type.y")] <- "Down_gene_type"
colnames(test)[which(colnames(test) == "gene_id.x")] <- "gene_id"
test <- test%>%mutate(Comparison = paste0(condition_1, " vs ", condition_2))
test <- test[, c(2, ncol(test), 6, 7, 10:(ncol(test)-1))]
write.csv(test, paste0(direOut,'consequence_detailed.csv'))


consequence_summary <- extractConsequenceSummary(
  asswitchlist_analyzed_final,
  consequencesToAnalyze=consequencesOfInterest,
  plotGenes = TRUE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE,      # enables analysis of fraction of significant features
  returnResult = TRUE
)

write.csv(consequence_summary,paste0(direOut,'consequence_summary.csv'))

consequence_comparison_switch <- extractConsequenceEnrichment(
  asswitchlist_analyzed_final,
  consequencesToAnalyze=consequencesOfInterest,
  analysisOppositeConsequence = TRUE,
  countGenes = FALSE,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
write.csv(consequence_comparison_switch,paste0(direOut,'consequence_comparison_switch.csv'))


consequence_comparison_gene <- extractConsequenceEnrichment(
  asswitchlist_analyzed_final,
  consequencesToAnalyze=consequencesOfInterest,
  analysisOppositeConsequence = TRUE,
  countGenes = TRUE,
  returnResult = TRUE # if TRUE returns a data.frame with the summary statistics
)
write.csv(consequence_comparison_gene,paste0(direOut,'consequence_comparison_gene.csv'))


##DEBUG: 5_utr_seq_similarity ####
consequencesOfInterest <- c('5_utr_seq_similarity')
tmp <- analyzeSwitchConsequences(asswitchlist_analyzed, consequencesToAnalyze=consequencesOfInterest,
                                 removeNonConseqSwitches = FALSE)

##DEBUG: analysis from scratch ####
pairwiseIsoComparison <- extractSwitchPairs(
  asswitchlist_analyzed,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = FALSE
)

pairwiseIsoComparisonUniq <-
  unique(pairwiseIsoComparison[,c(
    'isoformUpregulated', 'isoformDownregulated'
  )])
pairwiseIsoComparisonUniq$comparison <-
  1:nrow(pairwiseIsoComparisonUniq)


consequencesOfIsoformSwitching <- plyr::dlply(
  .data = pairwiseIsoComparisonUniq,
  .variables = 'comparison',
  .parallel = FALSE,
  .inform = TRUE,
  # .progress = progressBar,
  .fun = function(pairwiseIsoComparisonUniq) {
    # aDF <- pairwiseIsoComparisonUniq[246,]
    compareAnnotationOfTwoIsoforms(
      consequencesToAnalyze = consequencesOfInterest,
      upIso                 = pairwiseIsoComparisonUniq$isoformUpregulated,
      downIso               = pairwiseIsoComparisonUniq$isoformDownregulated,
      ntCutoff              = 50,
      ntFracCutoff          = NULL,
      ntJCsimCutoff         = 0.8,
      AaCutoff              = 10,
      AaFracCutoff          = 0.5,
      AaJCsimCutoff         = 0.9,
      addDescription        = TRUE,
      testInput             = FALSE # already done by this function
    )
  }
)

compareAnnotationOfTwoIsoforms <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = 'all',
    upIso,
    downIso,
    addDescription = TRUE,
    onlyRepportDifferent = FALSE,
    ntCutoff = 50,
    ntFracCutoff = NULL,
    ntJCsimCutoff = 0.8,
    AaCutoff = 10,
    AaFracCutoff = 0.5,
    AaJCsimCutoff = 0.9,
    testInput = TRUE
) {
  isoComparison <-
    data.frame(
      isoformUpregulated = upIso,
      isoformDownregulated = downIso,
      featureCompared = consequencesToAnalyze,
      isoformsDifferent = NA
    )
  if (addDescription) {
    isoComparison$switchConsequence <- NA
  }
  
  transcriptData <- asswitchlist_analyzed$orfAnalysis
  transcriptData <- transcriptData[which(transcriptData$isoform_id %in% c(upIso, downIso)), ]
  
  
  
  if (all(!is.na(transcriptData$orfTransciptStart))) {
    localUpNt   <-
      XVector::subseq(upNtSeq,
             1,
             transcriptData$orfTransciptStart[which(
               transcriptData$isoform_id == upIso
             )] - 1)
    localDownNt <-
      XVector::subseq(downNtSeq,
             1,
             transcriptData$orfTransciptStart[which(
               transcriptData$isoform_id == downIso
             )] -1)
    
    # if one of them have no UTR
    if (width(localUpNt) > 0 &
        width(localDownNt) > 0) {
      localAlignment <-
        Biostrings::pairwiseAlignment(
          pattern = localUpNt,
          subject = localDownNt,
          type = 'overlap'
        )
      
      overlapSize <- min(c(nchar(
        gsub(
          '-',
          '',
          as.character(localAlignment@subject)
        )
      ),
      nchar(
        gsub(
          '-',
          '',
          as.character(localAlignment@pattern)
        )
      )))
      totalWidth <-
        width(localAlignment@subject@unaligned) +
        width(localAlignment@pattern@unaligned) -
        overlapSize
      
      jcDist <- overlapSize / totalWidth
    } else {
      overlapSize <- 0
      totalWidth <-
        abs(width(localUpNt) - width(localDownNt))
      jcDist <- 0
    }
    
    differenttUTRoverlap <-
      jcDist < ntJCsimCutoff & totalWidth - overlapSize > ntCutoff
    
    # make repport
    localIndex <-
      which(
        isoComparison$featureCompared == '5_utr_seq_similarity'
      )
    isoComparison$isoformsDifferent[localIndex] <-
      differenttUTRoverlap
    
    # make description
    if (differenttUTRoverlap & addDescription) {
      utr5Gain <-
        transcriptData$isoform_id[which.max(
          transcriptData$orfTransciptStart
        )] == upIso
      
      if (utr5Gain) {
        isoComparison$switchConsequence[localIndex] <-
          '5UTR is longer'
      } else {
        isoComparison$switchConsequence[localIndex] <-
          '5UTR is shorter'
      }
      
    }
  }
}




##DEBUG: individual test case ####
upIso <- pairwiseIsoComparisonUniq$isoformUpregulated[374]
downIso <- pairwiseIsoComparisonUniq$isoformDownregulated[374]

isoformsToAnalyze <-  c(upIso, downIso)
names(isoformsToAnalyze) <- c('up', 'down')

ntSeq <-
  asswitchlist_analyzed$ntSequence[which(
    names(asswitchlist_analyzed$ntSequence) %in% c(upIso, downIso)
  )]
ntSeq

transcriptData <- asswitchlist_analyzed$orfAnalysis
transcriptData <- transcriptData[which(transcriptData$isoform_id %in% isoformsToAnalyze), ]



upNtSeq   <- ntSeq[upIso]
downNtSeq <- ntSeq[downIso]



localUpNt   <-
  XVector::subseq(upNtSeq,
         1,
         transcriptData$orfTransciptStart[which(
           transcriptData$isoform_id %in% upIso
         )] - 1)
localUpNt


localDownNt <-
  XVector::subseq(downNtSeq,
         1,
         transcriptData$orfTransciptStart[which(
           transcriptData$isoform_id %in% downIso
         )] -1)
localDownNt


# localAlignment <-
#   Biostrings::pairwiseAlignment(pattern = upNtSeq,
#                                 subject = downNtSeq,
#                                 type = 'global')
localAlignment <-
  Biostrings::pairwiseAlignment(
    pattern = localUpNt,
    subject = localDownNt,
    type = 'overlap'
  )
localAlignment@subject


# overlapSize <- min(c(nchar(gsub(
#   '-', '', as.character(localAlignment@subject)
# )),
# nchar(gsub(
#   '-', '', as.character(localAlignment@pattern)
# ))))
overlapSize <- min(c(nchar(
  gsub(
    '-',
    '',
    as.character(localAlignment@subject)
  )
),
nchar(
  gsub(
    '-',
    '',
    as.character(localAlignment@pattern)
  )
)))


overlapSize

totalWidth <-
  width(localAlignment@subject@unaligned) +
  width(localAlignment@pattern@unaligned) -
  overlapSize
totalWidth



###Alluvial plot ####
consequence_summary <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/consequence_summary.csv", row.names = 1)
consequence_summary$featureCompared <- gsub('\n', ' ', consequence_summary$featureCompared)

##Selected features ####
feature = c("NMD status", "ORF seq similarity", 'Domain length', 'Domains identified', 'Signal peptide identified')

tmp <- consequence_summary %>% filter(featureCompared %in% feature)
tmp_hyper <- tmp %>% group_by(Comparison, featureCompared) %>% summarise(nrGenes = sum(nrGenesWithConsequences))
# tmp <- tmp %>% group_by(Comparison, switchConsequence) %>% summarise(nrGenes = sum(nrGenesWithConsequences))
tmp <- tmp %>% filter(Comparison != "day0 vs day5")
tmp_hyper <- tmp_hyper %>% filter(Comparison != "day0 vs day5")

# is_lodes_form(tmp, key = "Comparison", value = "featureCompared", id = "nrGenes")

tmp <- tmp[order(tmp$featureCompared, tmp$nrGenesWithConsequences, decreasing = TRUE), ]
tmp_hyper <- tmp_hyper[order(tmp_hyper$Comparison, tmp_hyper$nrGenes, decreasing = TRUE), ]

level <- unique(tmp_hyper$featureCompared)
# level2 <- unique(tmp$switchConsequence)
level2 <- c("ORF is shorter", "ORF is longer", "Complete ORF loss", "Complete ORF gain", 
            'Domain loss', 'Domain gain', 'Domain switch',
            'NMD sensitive', 'NMD insensitive',
            'Domain length loss', 'Domain length gain', 
            'Signal peptide loss', 'Signal peptide gain')

factor(tmp$featureCompared, levels = level)


pdf(paste0(plot_dir_out, 'alluvial_consequence_summary_grouped.pdf'), 8, 6)
ggplot(data = tmp,
       # aes(x = Comparison, y = nrGenes, alluvium = featureCompared, stratum = featureCompared)) +
       aes(x = Comparison, y = nrGenesWithConsequences, alluvium = factor(switchConsequence, levels = level2),
           stratum = factor(featureCompared, levels = level))) +
       # aes(x = Comparison, y = nrGenesWithConsequences, alluvium = switchConsequence)) +
  theme_bw() + 
  geom_bar(aes(fill = factor(switchConsequence, levels = level2)),
           position = "stack")
  # scale_fill_gdocs()
  # scale_fill_brewer(type = "qual", palette = "gdocs")
  # proportional knot positioning (default)
  # geom_alluvium(aes(fill = factor(switchConsequence, levels = level2)),
  #               alpha = .75, decreasing = FALSE, width = 0.2) +
  # geom_stratum(aes(stratum = factor(featureCompared, levels = level)), decreasing = FALSE, width = 0.2) +
  # geom_alluvium(aes(fill = featureCompared),
  #               alpha = .75, decreasing = FALSE, width = 0.1) +
  # geom_stratum(aes(stratum = featureCompared), decreasing = FALSE, width = 0.1) +
  # scale_fill_gdocs()
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(title = "Functional Consequences", size = 10)) +
  xlab("State transition") +
  ylab("Number of genes")
dev.off()

pdf(paste0(plot_dir_out, 'alluvial_consequence_summary_ungrouped.pdf'), 8, 6)
ggplot(data = tmp,
       # aes(x = Comparison, y = nrGenes, alluvium = featureCompared, stratum = featureCompared)) +
       aes(x = Comparison, y = nrGenesWithConsequences, alluvium = factor(switchConsequence, levels = level2),
           stratum = factor(switchConsequence, levels = level2))) +
  # aes(x = Comparison, y = nrGenesWithConsequences, alluvium = switchConsequence)) +
  theme_bw() + 
  # scale_fill_gdocs()
  # scale_fill_brewer(type = "qual", palette = "gdocs")
  # proportional knot positioning (default)
  geom_alluvium(aes(fill = factor(switchConsequence, levels = level2)),
                alpha = .75, decreasing = FALSE, width = 0.2) +
  geom_stratum(aes(stratum = factor(switchConsequence, levels = level2)), decreasing = FALSE, width = 0.2) +
  # geom_alluvium(aes(fill = featureCompared),
  #               alpha = .75, decreasing = FALSE, width = 0.1) +
  # geom_stratum(aes(stratum = featureCompared), decreasing = FALSE, width = 0.1) +
  # scale_fill_gdocs()
  scale_fill_manual(values = colors) +
  guides(fill = guide_legend(title = "Functional Consequences", size = 10)) +
  xlab("State transition") +
  ylab("Number of genes")
dev.off()
##stacked bar plot ####
pdf(paste0(plot_dir_out, 'barplot_consequence_summary_grouped.pdf'), 5, 4)
ggplot(data = tmp,
       # aes(x = Comparison, y = nrGenes, alluvium = featureCompared, stratum = featureCompared)) +
       # aes(x = Comparison, y = nrGenesWithConsequences, alluvium = factor(switchConsequence, levels = level2),
       #     stratum = factor(featureCompared, levels = level))) +
      aes(x = Comparison, y = nrGenesWithConsequences)) +
  theme_bw() + 
  geom_bar(aes(fill = factor(switchConsequence, levels = level2)),
           position = "stack", stat="identity") + 
# scale_fill_gdocs()
# scale_fill_brewer(type = "qual", palette = "gdocs")
# proportional knot positioning (default)
# geom_alluvium(aes(fill = factor(switchConsequence, levels = level2)),
#               alpha = .75, decreasing = FALSE, width = 0.2) +
# geom_stratum(aes(stratum = factor(featureCompared, levels = level)), decreasing = FALSE, width = 0.2) +
# geom_alluvium(aes(fill = featureCompared),
#               alpha = .75, decreasing = FALSE, width = 0.1) +
# geom_stratum(aes(stratum = featureCompared), decreasing = FALSE, width = 0.1) +
# scale_fill_gdocs()
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = c("day0 vs day3" = "day 3 vs day 0",
                     "day3 vs day5" = "day 5 vs day 3")) + 
  guides(fill = guide_legend(title = "Functional Consequences", size = 10)) +
  # xlab("State transition") +
  ylab("Number of genes") +
  theme(text = element_text(size = 10, color = 'black'),
        axis.title.x = element_blank(), 
        axis.text = element_text(size = 10, color = 'black'))
dev.off()

##RNA-level ####
rna_feature = c("Tss", "Tts", "Exon number", "Intron structure", "Coding potential", "Isoform seq similarity", 
               "ORF length", "3 utr seq similarity", "5 utr length", "Last exon")
tmp <- consequence_summary %>% filter(featureCompared %in% rna_feature)
tmp <- tmp %>% group_by(Comparison, featureCompared) %>% summarise(nrGenes = sum(nrGenesWithConsequences))
tmp <- tmp %>% filter(Comparison != "day0 vs day5")
# is_lodes_form(tmp, key = "Comparison", value = "featureCompared", id = "nrGenes")

pdf(paste0(plot_dir_out, 'alluvial_consequence_summary_rna.pdf'), 8, 6)
ggplot(data = tmp,
             aes(x = Comparison, y = nrGenes, alluvium = featureCompared, stratum = featureCompared)) +
  theme_bw() + 
  # scale_fill_gdocs()
  # scale_fill_brewer(type = "qual", palette = "gdocs")
# proportional knot positioning (default)

  geom_alluvium(aes(fill = featureCompared),
                alpha = .75, decreasing = FALSE, width = 0.1) +
  geom_stratum(aes(stratum = featureCompared), decreasing = FALSE, width = 0.1) +
  # scale_fill_gdocs()
  scale_fill_manual(values = mycolors) +
  guides(fill = guide_legend(title = "RNA-level Functional Consequences", size = 10)) +
  xlab("State transition") +
  ylab("Number of genes")
dev.off()




##ORF-level ####
orf_feature = c("NMD status", "ORF seq similarity")
tmp <- consequence_summary %>% filter(featureCompared %in% orf_feature)
tmp <- tmp %>% group_by(Comparison, featureCompared) %>% summarise(nrGenes = sum(nrGenesWithConsequences))
tmp <- tmp %>% filter(Comparison != "day0 vs day5")
# is_lodes_form(tmp, key = "Comparison", value = "featureCompared", id = "nrGenes")
pdf(paste0(plot_dir_out, 'alluvial_consequence_summary_orf.pdf'), 8, 6)
ggplot(data = tmp,
       aes(x = Comparison, y = nrGenes, alluvium = featureCompared, stratum = featureCompared)) +
  theme_bw() + 
  # scale_fill_gdocs()
  # scale_fill_brewer(type = "qual", palette = "gdocs")
  # proportional knot positioning (default)
  
  geom_alluvium(aes(fill = featureCompared),
                alpha = .75, decreasing = FALSE, width = 0.1) +
  geom_stratum(aes(stratum = featureCompared), decreasing = FALSE, width = 0.1) +
  # scale_fill_gdocs()
  scale_fill_manual(values = mycolors) +
  guides(fill = guide_legend(title = "ORF-level Functional Consequences", size = 10)) +
  xlab("State transition") +
  ylab("Number of genes")
dev.off()

##Peptide-level ####
aa_feature = c('Domain length', 'Domains identified', 'Signal peptide identified')
tmp <- consequence_summary %>% filter(featureCompared %in% aa_feature)
tmp <- tmp %>% group_by(Comparison, featureCompared) %>% summarise(nrGenes = sum(nrGenesWithConsequences))
tmp <- tmp %>% filter(Comparison != "day0 vs day5")
# is_lodes_form(tmp, key = "Comparison", value = "featureCompared", id = "nrGenes")
pdf(paste0(plot_dir_out, 'alluvial_consequence_summary_aa.pdf'), 8, 6)
ggplot(data = tmp,
       aes(x = Comparison, y = nrGenes, alluvium = featureCompared, stratum = featureCompared)) +
  theme_bw() + 
  # scale_fill_gdocs()
  # scale_fill_brewer(type = "qual", palette = "gdocs")
  # proportional knot positioning (default)
  
  geom_alluvium(aes(fill = featureCompared),
                alpha = .75, decreasing = FALSE, width = 0.1) +
  geom_stratum(aes(stratum = featureCompared), decreasing = FALSE, width = 0.1) +
  # scale_fill_gdocs()
  scale_fill_manual(values = mycolors) +
  guides(fill = guide_legend(title = "Peptide-level Functional Consequences", size = 10)) +
  xlab("State transition") +
  ylab("Number of genes")
dev.off()




###isoform switch with multiple consequences ####
switch <- read.csv("/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/consequence_detailed.csv", row.names = 1)
switch <- switch %>% filter(Comparison != "day0 vs day5" & isoformsDifferent == TRUE)
test <- switch %>% filter(switchConsequence == "Domain loss")

tmp <- switch[switch$featureCompared == "domain_length", ]

NA %in% switch$gene_id

######
switch_summary<-as.data.frame(asswitchlist_analyzed$isoformFeatures)
write.csv(switch_summary,paste0(direOut,'switch_summary.csv'))

total_summary <- extractSwitchSummary(asswitchlist_analyzed)
write.csv(total_summary,paste0(direOut,'total_summary.csv'))


switch <- extractTopSwitches(
  extractGenes=TRUE, 
  asswitchlist_analyzed, 
  n = NA, 
  sortByQvals = TRUE
)
write.csv(switch,paste0(direOut,'switch_gene_list.csv'))




asswitchlist_analyzed <- analyzeAlternativeSplicing(switchAnalyzeRlist = asswitchlist_analyzed)

event_summary<-as.data.frame(asswitchlist_analyzed$AlternativeSplicingAnalysis)
# write.csv(event_summary,paste0(direOut,'AS_type_detail.csv'))

df <- extractSplicingSummary(asswitchlist_AS_analyzed, asFractionTotal = FALSE, plotGenes=FALSE, plot = FALSE, returnResult = TRUE)
# write.csv(df,paste0(direOut,'AS_type_summary.csv'))

splicing_comparison <- extractSplicingEnrichment(
  asswitchlist_analyzed,
  splicingToAnalyze='all',
  countGenes = FALSE,
  returnResult=TRUE,
  returnSummary=TRUE
)
write.csv(splicing_comparison,paste0(direOut,'splicing_comparison.csv'))

splicing_summary <- extractSplicingGenomeWide(asswitchlist_analyzed,plot=TRUE)
write.csv(splicing_summary,paste0(direOut,'splicing_summary.csv'))











top <- extractTopSwitches(asswitchlist_analyzed, 
                   filterForConsequences = FALSE, 
                   sortByQvals = TRUE,
                   n = NA)
row.names(top) <- top$Rank
###############
#????????
###############
# The next two lines were done to filter out misannotated genes
# htop <- table(top$gene_switch_q_value)
# top <- top[top$gene_switch_q_value %in% names(htop[which(htop == 1)]),] 
###############



sig_genes <- unique(top$gene_name)



extractSwitchOverlap(
  asswitchlist_analyzed,
  filterForConsequences=FALSE,scaleVennIfPossible = FALSE
)

isoPairs <- extractSwitchPairs(
  switchAnalyzeRlist = asswitchlist_analyzed,
  alpha = 0.05,
  dIFcutoff = 0.1,
  onlySigIsoforms = TRUE
)
isoPairs$switch <- stringr::str_c(
  isoPairs$isoformDownregulated,
  isoPairs$isoformUpregulated
)

isoPairs$comparison <- stringr::str_c(
  isoPairs$condition_1,
  '\nvs\n',
  isoPairs$condition_2
)
switchList <- split(isoPairs$switch, isoPairs$comparison)





switchPlot(asswitchlist_analyzed, gene = 'PFN2', condition1 = "day0", condition2 = "day5")

splicingEnrichment <- extractSplicingEnrichment(
  asswitchlist_analyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)


extractSwitchOverlap(asswitchlist_analyzed)



extractSwitchPairs <- function(
    switchAnalyzeRlist,
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE
) {
  ### Check input
  if (TRUE) {
    # check switchAnalyzeRlist
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
      stop(
        'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
      )
    }
    
    if (alpha < 0 |
        alpha > 1) {
      warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
    }
    if (alpha > 0.05) {
      warning(
        'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
      )
    }
    
    # test wether switching have been analyzed
    if (!any(!is.na(
      switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
    ))) {
      stop(
        'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.'
      )
    }
  }
  
  ### Extract and massage data
  if (TRUE) {
    localData <- switchAnalyzeRlist$isoformFeatures[
      which(
        switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
          abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
      ),
      c(
        'iso_ref',
        'gene_ref',
        'isoform_switch_q_value',
        'gene_switch_q_value',
        'dIF'
      )
    ]
    
    if (!nrow(localData)) {
      stop('No genes were considered switching with the used cutoff values')
    }
    
    ### add switch direction
    localData$switchDirection <- NA
    localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
    localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'
    
    ### Annotate significant features
    isoResTest <-
      any(!is.na(
        switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
      ))
    if (isoResTest) {
      localData$isoSig <-
        localData$isoform_switch_q_value < alpha &
        abs(localData$dIF) > dIFcutoff
    } else {
      localData$isoSig <-
        localData$gene_switch_q_value < alpha &
        abs(localData$dIF) > dIFcutoff
    }
    
    
    if(onlySigIsoforms) {
      localData <- localData[which( localData$isoSig ),]
    }
  }
  
  
  ### Create data sub-sets of interest
  if(TRUE) {
    sigUpData <- localData[which(
      localData$isoSig & localData$switchDirection == 'up'
    ),c('iso_ref','gene_ref')]
    sigDnData <- localData[which(
      localData$isoSig & localData$switchDirection == 'down'
    ),c('iso_ref','gene_ref')]
    
    colnames(sigUpData)[1] <- c('iso_ref_up')
    colnames(sigDnData)[1] <- c('iso_ref_down')
    
    
    if( ! onlySigIsoforms ) {
      justUpData <- localData[which(
        localData$switchDirection == 'up'
      ),c('iso_ref','gene_ref')]
      justDnData <- localData[which(
        localData$switchDirection == 'down'
      ),c('iso_ref','gene_ref')]
      
      colnames(justUpData)[1]  <- c('iso_ref_up')
      colnames(justDnData)[1] <- c('iso_ref_down')
    }
  }
  
  ### Join datasets to extract pairs
  if(TRUE) {
    if( onlySigIsoforms ) {
      pairwiseIsoComparison <- dplyr::inner_join(
        sigUpData,
        sigDnData,
        by= 'gene_ref'
      )
    } else {
      ### Sig up and all down
      upPairs <- dplyr::inner_join(
        sigUpData,
        justDnData,
        by= 'gene_ref'
      )
      ### Sig down and all up
      dnPairs <- dplyr::inner_join(
        justUpData,
        sigDnData,
        by= 'gene_ref'
      )
      
      ### Combine
      pairwiseIsoComparison <- unique(
        rbind(
          upPairs,
          dnPairs
        )
      )
      pairwiseIsoComparison <- pairwiseIsoComparison[,c(
        'gene_ref','iso_ref_up','iso_ref_down'
      )]
      
      ### Reorder
      pairwiseIsoComparison <- pairwiseIsoComparison[order(
        pairwiseIsoComparison$gene_ref,
        pairwiseIsoComparison$iso_ref_up,
        pairwiseIsoComparison$iso_ref_down
      ),]
    }
  }
  
  ### Add in additional data
  if(TRUE) {
    ### Add isoform names
    matchVectorUp <- match(
      pairwiseIsoComparison$iso_ref_up,
      switchAnalyzeRlist$isoformFeatures$iso_ref
    )
    matchVectorDn <- match(
      pairwiseIsoComparison$iso_ref_down,
      switchAnalyzeRlist$isoformFeatures$iso_ref
    )
    
    pairwiseIsoComparison$isoformUpregulated   <-
      switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorUp]
    pairwiseIsoComparison$isoformDownregulated <-
      switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorDn]
    
    ### Gene infl
    pairwiseIsoComparison$gene_id <-
      switchAnalyzeRlist$isoformFeatures$gene_id[matchVectorUp]
    pairwiseIsoComparison$gene_name <-
      switchAnalyzeRlist$isoformFeatures$gene_name[matchVectorUp]
    
    ### Conditons
    pairwiseIsoComparison$condition_1 <-
      switchAnalyzeRlist$isoformFeatures$condition_1[matchVectorUp]
    pairwiseIsoComparison$condition_2 <-
      switchAnalyzeRlist$isoformFeatures$condition_2[matchVectorUp]
    
  }
  
  return( pairwiseIsoComparison )
}


# splicingToAnalyze = 'all'
# asFractionTotal = FALSE
# alpha = 0.05
# dIFcutoff = 0.1
# onlySigIsoforms = FALSE
# plot = TRUE
# plotGenes = FALSE
# localTheme = theme_bw()
# returnResult = FALSE
# 
# switchAnalyzeRlist <- asswitchlist_analyzed
# localData <- switchAnalyzeRlist$isoformFeatures[which(
#   switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < alpha &
#     abs(switchAnalyzeRlist$isoformFeatures$dIF) > dIFcutoff
# ),
# c(
#   'iso_ref',
#   'gene_ref',
#   'isoform_switch_q_value',
#   'gene_switch_q_value',
#   'dIF'
# )]
# localData$switchDirection <- NA
# localData$switchDirection[which(sign(localData$dIF) ==  1)] <- 'up'
# localData$switchDirection[which(sign(localData$dIF) == -1)] <- 'down'
# 
# ### split based on genes and conditions
# localDataList <-
#   split(localData, f = localData$gene_ref, drop = TRUE)
# 
# 
# pairwiseIsoComparisonList <-
#   plyr::llply(
#     .data = localDataList,
#     .progress = 'none',
#     .fun = function(aDF) {
#       # aDF <- localDataList[[171]]
#       isoResTest <-
#         any(!is.na(
#           switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value
#         ))
#       if (isoResTest) {
#         sigIso <- aDF$iso_ref[which(
#           aDF$isoform_switch_q_value < alpha &
#             abs(aDF$dIF) > dIFcutoff
#         )]
#       } else {
#         sigIso <- aDF$iso_ref[which(
#           aDF$gene_switch_q_value < alpha &
#             abs(aDF$dIF) > dIFcutoff
#         )]
#       }
#       if (length(sigIso) == 0) {
#         return(NULL)
#       }
#       
#       ### reduce to significant if nessesary
#       if (onlySigIsoforms) {
#         aDF <- aDF[which(aDF$iso_ref %in% sigIso), ]
#       }
#       if (nrow(aDF) < 2) {
#         return(NULL)
#       }
#       
#       ### make sure there are both up and down
#       if (!all(c('up', 'down') %in% aDF$switchDirection)) {
#         return(NULL)
#       }
#       
#       ### extract pairs of isoforms
#       upIso   <-
#         as.vector(aDF$iso_ref[which(
#           aDF$switchDirection == 'up'
#         )])
#       downIso <-
#         as.vector(aDF$iso_ref[which(
#           aDF$switchDirection == 'down'
#         )])
#       
#       allIsoCombinations <-
#         setNames(
#           base::expand.grid(
#             upIso,
#             downIso,
#             stringsAsFactors = FALSE,
#             KEEP.OUT.ATTRS = FALSE
#           ),
#           nm = c('iso_ref_up', 'iso_ref_down')
#         )
#       
#       ### Reduce to those where at least one of them is significant
#       allIsoCombinations <-
#         allIsoCombinations[which(
#           allIsoCombinations$iso_ref_up %in% sigIso |
#             allIsoCombinations$iso_ref_down %in% sigIso
#         ), ]
#       
#       ### Add gen ref
#       allIsoCombinations$gene_ref    <- aDF$gene_ref[1]
#       
#       return(allIsoCombinations)
#     }
#   )
# 
# pairwiseIsoComparisonList <-
#   pairwiseIsoComparisonList[which(
#     ! sapply(pairwiseIsoComparisonList, is.null)
#   )]
# 
# myListToDf <- function(
#     aList, # List with data.frames to concatenate
#     ignoreColNames = FALSE, # A Logical indicating whether to check the colnames of each data.frame in aList
#     addOrignAsRowNames = FALSE, # A Logical indicating whether to add the name of the list intry as rownames in the final data.frame
#     addOrignAsColumn = FALSE, # A logical indicating whether a column conatining the name of the list entry should be added in the final data.frame
#     addOrgRownames = FALSE # A logical indicating whther the original rownames should be used in the final data.frame
# ) {
#   ### Test whether input match standards for being bound together
#   if (class(aList) != 'list') {
#     stop("Input is not a list")
#   }
#   
#   # remove empty ones
#   aList <- aList[which(!sapply(aList, is.null))]
#   
#   # Make sure the list entries are data.frames
#   if (class(aList[[1]]) != "data.frame") {
#     aList <- lapply(aList, function(x)
#       as.data.frame(t(x)))
#   }
#   
#   nCol <- unique(sapply(aList, ncol))
#   if (length(nCol)  != 1) {
#     stop("Interies in the list does not have the same number of collums/")
#   }
#   if (!ignoreColNames) {
#     if (length(unique(as.vector(sapply(
#       aList, names
#     )))) !=  nCol) {
#       stop("Interies in the list does not have the collum names")
#     }
#   }
#   
#   ### data.frame to store results
#   df <-
#     data.frame(matrix(NA, ncol = nCol, nrow = sum(sapply(aList, nrow))))
#   
#   ### use sapply to loop over the list and extract the entries one at the time
#   for (i in 1:nCol) {
#     df[, i] <-
#       as.vector(unlist(sapply(aList, function(x)
#         x[, i]))) # the combination of as.vector and unlist makes it posible to have any number of entries in each of the lists
#   }
#   
#   # add names
#   colnames(df) <- colnames(aList[[1]])
#   if (addOrignAsColumn)    {
#     df$orign     <- rep(names(aList)            , sapply(aList, nrow))
#   }
#   if (addOrignAsRowNames)  {
#     rownames(df) <- rep(names(aList)            , sapply(aList, nrow))
#   }
#   if (addOrgRownames)      {
#     rownames(df) <- rep(sapply(aList, rownames) , sapply(aList, nrow))
#   }
#   
#   return(df)
# }
# pairwiseIsoComparison <-
#   myListToDf(pairwiseIsoComparisonList, addOrignAsColumn = FALSE)
# 
# pairwiseIsoComparison$isoformUpregulated   <-
#   switchAnalyzeRlist$isoformFeatures$isoform_id[match(
#     pairwiseIsoComparison$iso_ref_up,
#     switchAnalyzeRlist$isoformFeatures$iso_ref
#   )]
# pairwiseIsoComparison$isoformDownregulated <-
#   switchAnalyzeRlist$isoformFeatures$isoform_id[match(
#     pairwiseIsoComparison$iso_ref_down,
#     switchAnalyzeRlist$isoformFeatures$iso_ref
#   )]
# 
# # gene info
# pairwiseIsoComparison$gene_id   <-
#   switchAnalyzeRlist$isoformFeatures$gene_id[match(
#     pairwiseIsoComparison$iso_ref_up,
#     switchAnalyzeRlist$isoformFeatures$iso_ref
#   )]
# pairwiseIsoComparison$gene_name   <-
#   switchAnalyzeRlist$isoformFeatures$gene_name[match(
#     pairwiseIsoComparison$iso_ref_up,
#     switchAnalyzeRlist$isoformFeatures$iso_ref
#   )]
# # condition
# pairwiseIsoComparison$condition_1 <-
#   switchAnalyzeRlist$isoformFeatures$condition_1[match(
#     pairwiseIsoComparison$iso_ref_down,
#     switchAnalyzeRlist$isoformFeatures$iso_ref
#   )]
# 
# pairwiseIsoComparison$condition_2 <-
#   switchAnalyzeRlist$isoformFeatures$condition_2[match(
#     pairwiseIsoComparison$iso_ref_down,
#     switchAnalyzeRlist$isoformFeatures$iso_ref
#   )]
# 
# localAS <- switchAnalyzeRlist$AlternativeSplicingAnalysis
# localAS <- localAS[which(
#   localAS$isoform_id %in% pairwiseIsoComparison$isoformUpregulated |
#     localAS$isoform_id %in% pairwiseIsoComparison$isoformDownregulated
# ),]
# 
# ### Massage
# m1 <- reshape2::melt(localAS[,c(
#   "isoform_id",
#   "ES_genomic_start",
#   "MEE_genomic_start",
#   "MES_genomic_start",
#   "IR_genomic_start",
#   "A5_genomic_start",
#   "A3_genomic_start",
#   "ATSS_genomic_start",
#   "ATTS_genomic_start"
# )], id.vars = 'isoform_id')
# colnames(m1)[3] <- 'genomic_start'
# m1$AStype <- sapply(
#   strsplit(as.character(m1$variable),'_'),
#   function(x) x[1]
# )
# 
# ### Add AS to pairs
# localConseq2 <- merge(
#   pairwiseIsoComparison,
#   m1[,c('isoform_id','AStype','genomic_start')],
#   by.x='isoformUpregulated',
#   by.y='isoform_id'
# )
# 
# localConseq3 <- merge(
#   localConseq2,
#   m1[,c('isoform_id','AStype','genomic_start')],
#   by.x=c('isoformDownregulated','AStype'),
#   by.y=c('isoform_id','AStype'),
#   suffixes = c("_up","_down")
# )
# 
# ### Summarize
# localConseq3$upAs <- !is.na(localConseq3$genomic_start_up)
# localConseq3$dnAs <- !is.na(localConseq3$genomic_start_down)
# 
# localConseq3$Comparison <-
#   paste(
#     localConseq3$condition_1,
#     'vs',
#     localConseq3$condition_2,
#     sep = ' '
#   )
# 
# tmpUp <- localConseq3[c('isoformUpregulated','gene_id','Comparison','AStype','upAs')]
# tmpDn <- localConseq3[c('isoformDownregulated','gene_id','Comparison','AStype','dnAs')]
# 
# colnames(tmpUp) <- c('isoform_id','gene_id','Comparison','AStype','anyAs')
# colnames(tmpDn) <- c('isoform_id','gene_id','Comparison','AStype','anyAs')
# 
# tmpUp$isoRegulation <- 'Up'
# tmpDn$isoRegulation <- 'Dn'
# 
# localConseq5 <- rbind(
#   tmpUp,
#   tmpDn
# )
# 
# myNumbers <- plyr::ddply(
#   localConseq5,
#   .variables = c(
#     'AStype',
#     'isoRegulation',
#     'Comparison'
#   ),
#   .drop = TRUE,
#   .fun = function(
#     aDF
#   ) { # aDF <- localConseq5[which( localConseq5$AStype == 'IR' & localConseq5$Comparison == 'COAD_ctrl vs COAD_cancer' & localConseq5$isoRegulation == 'Up'),]
#     localRes <- data.frame(
#       nrGenesWithConsequences = length(
#         unique( aDF$gene_id[which(aDF$anyAs)] )
#       ),
#       nrIsoWithConsequences = length(
#         unique( aDF$isoform_id[which(aDF$anyAs)] )
#       ),
#       stringsAsFactors = FALSE
#     )
#     
#     if (asFractionTotal) {
#       localRes$geneFraction <- localRes$nrGenesWithConsequences /
#         length( unique( aDF$gene_id ) )
#       localRes$isoFraction <- localRes$nrIsoWithConsequences /
#         length( unique( aDF$isoform_id ) )
#     }
#     return(localRes)
#   }
# )
# 
# myNumbers$splicingResult <- paste(
#   myNumbers$AStype,
#   ifelse(myNumbers$isoRegulation == 'Dn','in isoform used less','in isoform used more'),
#   sep=' '
# )
# 
# 
# 
# write.csv(pairwiseIsoComparison, paste0(direOut,'paired_isoswitch.csv'))
# write.csv(localConseq5,paste0(direOut,'isoform_AS_byDay.csv'))
# 

