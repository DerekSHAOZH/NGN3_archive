list.of.packages <- c("ggplot2","rjson")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("rjson"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("stringr"))

args <- commandArgs(trailingOnly = TRUE) 

outDir = args[1]
dataDir = args[2]

getFigure = function(GO.path, FDR_threshold, breaks){
  level = 0 # <----- change plot level if desired 
  enrichment_type="+" # <----- change enrichment type to "-" if you want to plot underrepresentet GO terms 
  
  result <- fromJSON(file = GO.path)
  
  
  shortCut = result$overrepresentation$group
  TABLE=data.frame()
  j=0
  for (level0 in 1:length(shortCut)){
    i=0
    j=j+1
    if(length(shortCut)==1){
      shortCut_result=list()
      shortCut_result[[1]]=shortCut$result
      
    }else{
      if(is.null(shortCut[[level0]]$result[[1]]$term$level)){
        shortCut_result=list()
        shortCut_result[[1]]=shortCut[[level0]]$result
      }else{
        shortCut_result=shortCut[[level0]]$result
      }
    }
    sublevel=length(shortCut_result)
    notIncluded=TRUE
    for (level_other in 1:sublevel){
      
      level_obj=NULL
      label_obj=NULL
      ROW=NULL
      
      
      if (!is.null(shortCut_result[[level_other]]$term$level) && !is.null(shortCut_result[[level_other]]$term$label)){
        
        level_obj=shortCut_result[[level_other]]$term$level
        label_obj=shortCut_result[[level_other]]$term$label
        ROW=shortCut_result[[level_other]]$input_list
        
      }else{
        level_obj=0
        label_obj=shortCut_result[[level_other]]$term$label
        ROW=shortCut_result[[level_other]]$input_list
      }
      
      if ((level_obj<=level && ROW$plus_minus==enrichment_type && as.numeric(ROW$fdr)<=FDR_threshold) || 
          (level_obj>level && ROW$plus_minus==enrichment_type && as.numeric(ROW$fdr)<=FDR_threshold && notIncluded)){
        
        if(notIncluded){
          name=label_obj
        }else{
          name=paste(rep(" -> ",level_obj),label_obj,sep="")
        }
        notIncluded=FALSE
        i=i+1
        
        observed=as.numeric(ROW$number_in_list)
        enrichment=as.numeric(ROW$fold_enrichment)
        FDR=as.numeric(ROW$fdr)
        intra_group_order=i
        extra_group_order=j
        
        TABLE=rbind(TABLE,
                    data.frame(name=name,observed=observed,
                               enrichment=enrichment,
                               FDR=FDR,
                               intra_group_order=intra_group_order,
                               extra_group_order=extra_group_order))
      }
    }
  }
  
  TABLE=unique(TABLE)
  
  new_order=order(TABLE$FDR[TABLE$intra_group_order==1])
  TABLE$extra_group_order_new=rep(pmatch(1:length(unique(TABLE$extra_group_order)),new_order),as.vector(table(TABLE$extra_group_order)))
  
  TABLE$name=gsub("via.*","via ...",TABLE$name)
  TABLE$name=gsub(",.*",", ...",TABLE$name)
  TABLE$name=gsub("positive ","pos. ",TABLE$name)
  TABLE$name=gsub("negative ","neg. ",TABLE$name)
  TABLE$name=stringr::str_wrap(TABLE$name,28)
  maxEnrichment = max(TABLE$enrichment)
  print(maxEnrichment)
  p = ggplot(data = subset(TABLE, FDR<FDR_threshold),aes(x=FDR, 
                                                     y = reorder(name,-1*extra_group_order_new),
                                                     #y=reorder(name,-1*FDR), 
                                                     colour=observed, 
                                                     size=enrichment)) +
    geom_point() +
    scale_size_continuous(range  = c(5, 20), 
                          limits = c(1,maxEnrichment), 
                          breaks = breaks) + 
    labs(x="FDR", y="", colour="DE genes", size="Enrichment") +
    theme_bw() +
    xlim(0, FDR_threshold) + 
    theme(#axis.text.y = element_text(hjust=0),
          axis.text=element_text(size=24),
          axis.title=element_text(size=24),
          legend.title=element_text(size=24),
          legend.text=element_text(size=24),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks.y=element_blank()
    )
  #dev.off()
  return(p)
}

FDR_threshold=0.01
breaks=c(1, 5, 8)

p1 = getFigure(paste(dataDir,"BP.day1vsday0.logFC1.up.json",sep="/"),FDR_threshold,breaks)
p2 = getFigure(paste(dataDir,"BP.day2vsday1.logFC1.up.json",sep="/"),FDR_threshold,breaks)
p3 = getFigure(paste(dataDir,"BP.day3vsday2.logFC1.up.json",sep="/"),FDR_threshold,breaks)
p4 = getFigure(paste(dataDir,"BP.day5vsday3.logFC1.up.json",sep="/"),FDR_threshold,breaks)


pdf(paste(outDir,"Figure_EV2B.pdf",sep="/"), width=35, height=20) #13
ggarrange(p1,p2,p3,p4,
          nrow = 1, ncol=4,
          align="hv", widths=c(1.5,1.5,1.5,1.5),
          common.legend = TRUE, legend="top")
dev.off()


p1 = getFigure(paste(dataDir,"BP.day1vsday0.logFC1.down.json",sep="/"),FDR_threshold,breaks)
p2 = getFigure(paste(dataDir,"BP.day2vsday1.logFC1.down.json",sep="/"),FDR_threshold,breaks)
p3 = getFigure(paste(dataDir,"BP.day3vsday2.logFC1.down.json",sep="/"),FDR_threshold,breaks)
p4 = getFigure(paste(dataDir,"BP.day5vsday3.logFC1.down.json",sep="/"),FDR_threshold,breaks)


pdf(paste(outDir,"Figure_EV2C.pdf",sep="/"), width=35, height=20) #13
ggarrange(p1,p2,p3,p4,
          nrow = 1, ncol=4,
          align="hv", widths=c(1.5,1.5,1.5,1.5),
          common.legend = TRUE, legend="top")
dev.off()

FDR_threshold=0.05
breaks=c(1,2,4)
p1 = getFigure(paste(dataDir,"Isoform.day1vsday0.fullGO.genericBackground.json",sep="/"),FDR_threshold,breaks)
p2 = getFigure(paste(dataDir,"Isoform.day3vsday2.fullGO.genericBackground.json",sep="/"),FDR_threshold,breaks)
p3 = getFigure(paste(dataDir,"Isoform.day5vsday3.fullGO.genericBackground.json",sep="/"),FDR_threshold,breaks)

pdf(paste(outDir,"Figure_3G.pdf",sep="/"), width=35, height=30) #13
ggarrange(p1,p2,p3,
          nrow = 1, ncol=3,
          align="hv", 
          common.legend = TRUE, legend="top")
dev.off()

