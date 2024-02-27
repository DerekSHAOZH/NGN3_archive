#BiocManager::install("IsoformSwitchAnalyzeR")
library("IsoformSwitchAnalyzeR")
library(ggplot2)
library(dplyr)

#set seed for reproductive results
set.seed(1) 

annotfn <- '/project/owlmayer/Applications/06_IsoformPipeline/data/genomeGencode/gencode.v28.annotation.gtf'
fastan <- '/project/owlmayer/Applications/06_IsoformPipeline/data/genomeGencode/gencode.transcripts.fa'
dire<-'/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Results/'
direout<-'/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitchAnalyser/Results/'
res_dir_out<-'/project/Neurodifferentiation_System/owlmayerTemporary/derek/isoSwitch/Results/'
samples <- c('OJ63','OJ64','OJ65','OJ66','OJ67','OJ68','OJ69','OJ70','OJ71','OJ72','OJ73','OJ74',
             'OJ75','OJ76','OJ77','OJ79','OJ80','OJ81','OJ83','OJ85','OJ86','OJ87','OJ88','OJ89')
cond <- c('day0','day0','day0','day0','day0','day1','day1','day1','day1','day1','day2','day2',
          'day2','day2','day2','day3','day3','day3','day3','day5','day5','day5','day5','day5')
rep <- c('rep1','rep2','rep3','rep4','rep5','rep1','rep2','rep3','rep4','rep5','rep1','rep2',
         'rep3','rep4','rep5','rep1','rep2','rep3','rep5','rep1','rep2','rep4','rep5','rep6')
fns <- paste0(dire,samples,'.isoforms.results')
names(fns)<-paste(cond,rep,sep='_')
myDesign <- data.frame(sampleID = paste(cond,rep,sep='_'),condition = cond)

rsemQuant <- importIsoformExpression(sampleVector=fns)



asswitchlist <- importRdata(
  isoformCountMatrix   = rsemQuant$counts,
  isoformRepExpression = rsemQuant$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = annotfn,
  isoformNtFasta       = fastan,
  showProgress = FALSE
)
tmp <- asswitchlist$isoformFeatures
# switchList <- isoformSwitchAnalysisPart1(
#   switchAnalyzeRlist   = asswitchlist,
#   pathToOutput = direout,
#   outputSequences      = TRUE,  
#   prepareForWebServers = FALSE)
#   

asswitchlist_filtered <- preFilter(switchAnalyzeRlist = asswitchlist,geneExpressionCutoff = 1,isoformExpressionCutoff = 0,removeSingleIsoformGenes = TRUE)
# test <- asswitchlist_filtered$isoformFeatures
# filtered_gene <- unique(test$gene_id)
# write.csv(filtered_gene,paste0(res_dir_out,'isoformswitchanalyzer_prefilter_gene.csv'), quote = F)

# asswitchlist_analyzed <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = asswitchlist_filtered, reduceToSwitchingGenes=TRUE)
asswitchlist_analyzed <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = asswitchlist_filtered, reduceToSwitchingGenes=FALSE)

extractSwitchSummary(asswitchlist_analyzed)
extractSwitchSummary(asswitchlist_analyzed, onlySigIsoforms = T)
switch <- extractTopSwitches(
  extractGenes=FALSE, 
  asswitchlist_analyzed, 
  n = NA, 
  sortByQvals = TRUE
)





switch_summary<-as.data.frame(asswitchlist_analyzed$isoformFeatures)
# write.csv(switch_summary,paste0(res_dir_out,'switch_summary.csv'))
write.csv(switch_summary,paste0(res_dir_out,'switch_summary.full.csv'))


test <- analyzeAlternativeSplicing(asswitchlist_analyzed)
event_summary<-as.data.frame(test$AlternativeSplicingAnalysis)
write.csv(event_summary,paste0(res_dir_out,'AS_type_detail.csv'))



df <- extractSplicingSummary(test, asFractionTotal = FALSE, plotGenes=FALSE, plot = FALSE, returnResult = TRUE)
write.csv(df,paste0(res_dir_out,'AS_type_summary.csv'))

asswitchlist_analyzed <- analyzeAlternativeSplicing(asswitchlist_analyzed)
asswitchlist_analyzed <- extractSequence(asswitchlist_analyzed)
asswitchlist_analyzed <- analyzeSwitchConsequences(asswitchlist_analyzed, consequencesToAnalyze=c(
  'intron_retention',
  'ORF_seq_similarity',
  'ORF_length',
  'NMD_status'
))


switchPlot(asswitchlist_analyzed, gene = 'TSC22D1', condition1 = "day1", condition2 = "day2")

switchPlotTopSwitches(
  switchAnalyzeRlist = asswitchlist_analyzed, 
  n = Inf,                                             # Set to Inf for all
  filterForConsequences = FALSE,
  fileType = "pdf",                                   # alternative is "png"
  pathToOutput = plot_dir_out
)

extractSplicingEnrichment(
  asswitchlist_analyzed,
  splicingToAnalyze='all',
  returnResult=FALSE,
  returnSummary=TRUE
)

extractConsequenceSummary(
  asswitchlist_analyzed,
  consequencesToAnalyze=c(
    'intron_retention',
    'ORF_seq_similarity',
    'ORF_length',
    'NMD_status'
  ),
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)

extractConsequenceEnrichment(
  asswitchlist_analyzed,
  consequencesToAnalyze=c(
    'intron_retention',
    'ORF_seq_similarity',
    'ORF_length',
    'NMD_status'
  ),
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)


##########################################################################
### only subsequent pairwise comparison
##########################################################################
days <- c(0, 1, 2, 3, 5)

for (i in 2:length(days)-1) {
  index <- which(grepl(days[i], cond, fixed = T))
  index <- c(index, which(grepl(days[i+1], cond, fixed = T)))
  print(index)
  
  samples_tmp <- samples[index]
  cond_tmp <- cond[index]
  rep_tmp <- rep[index]
  fns <- paste0(dire,samples_tmp,'.isoforms.results')
  names(fns)<-paste(cond_tmp,rep_tmp,sep='_')
  myDesign <- data.frame(sampleID = paste(cond_tmp,rep_tmp,sep='_'),condition = cond_tmp)
  
  rsemQuant <- importIsoformExpression(sampleVector=fns)
  
  asswitchlist <- importRdata(
    isoformCountMatrix   = rsemQuant$counts,
    isoformRepExpression = rsemQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = annotfn,
    isoformNtFasta       = fastan,
    showProgress = FALSE
  )
  
  asswitchlist <- preFilter(switchAnalyzeRlist = asswitchlist,geneExpressionCutoff = 1,isoformExpressionCutoff = 0,removeSingleIsoformGenes = TRUE)
  asswitchlist_analyzed <- isoformSwitchTestDEXSeq(switchAnalyzeRlist = asswitchlist, reduceToSwitchingGenes=TRUE)
  
  switch_summary<-as.data.frame(asswitchlist_analyzed$isoformFeatures)
  write.csv(switch_summary,paste0(res_dir_out,paste0('switch_summary_day', days[i], '_vs_day', days[i+1], '.csv')))
  
  
  test <- analyzeAlternativeSplicing(asswitchlist_analyzed)
  event_summary<-as.data.frame(test$AlternativeSplicingAnalysis)
  write.csv(event_summary,paste0(res_dir_out,paste0('AS_type_detail_day', days[i], '_vs_day', days[i+1], '.csv')))
  
  
  
  df <- extractSplicingSummary(test, asFractionTotal = FALSE, plotGenes=FALSE, plot = FALSE, returnResult = TRUE)
  write.csv(df,paste0(res_dir_out,paste0('AS_type_summary_day',days[i], '_vs_day', days[i+1], '.csv')))
}

