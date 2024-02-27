list.of.packages <- c("ggplot2","rjson")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("rjson"))

GO.path = "/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/RBPs/cluster_4.GOslim_Biological_Process.json" # <----- change jason file input
level = 0 # <----- change plot level if desired 
height =10 # <----- change plot heigth for pdf output
enrichment_type="+" # <----- change enrichment type to "-" if you want to plot underrepresentet GO terms 
out = "/project/owlmayer/Annkatrin/Share/GO/test.pdf" # <----- change output file

args <- commandArgs(trailingOnly = TRUE) 


GO.path = args[1]
level = as.numeric(args[2]) # <----- change plot level if desired 
height = as.numeric(args[3]) # <----- change plot heigth for pdf output
enrichment_type=args[4] # <----- change enrichment type to "-" if you want to plot underrepresentet GO terms 
out = args[5] # <----- change output file


result <- fromJSON(file = GO.path)


shortCut = result$overrepresentation$group
TABLE=data.frame()
j=0
for (level0 in 1:length(shortCut)){
  i=0
  j=j+1
  for (level_other in 1:length(shortCut[[level0]]$result)){
    
    if (is.null(shortCut[[level0]]$result$term$label)){
      level_obj = shortCut[[level0]]$result[[level_other]]$term$level
      label_obj = shortCut[[level0]]$result[[level_other]]$term$label
      ROW=shortCut[[level0]]$result[[level_other]]$input_list
      
    }else{
      level_obj = ifelse(is.null(shortCut[[level0]]$result$term$level),0,shortCut[[level0]]$result$term$level)
      label_obj = shortCut[[level0]]$result$term$label
      
      ROW=shortCut[[level0]]$result$input_list
      i=0
    }

    if (level_obj<=level && 
        ROW$plus_minus==enrichment_type){
      i=i+1
      name=paste(rep(" -> ",level_obj),label_obj,sep="")
      
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
TABLE$extra_group_order_new=rep(pmatch(1:max(TABLE$extra_group_order),new_order),as.vector(table(TABLE$extra_group_order)))

TABLE
max=max(TABLE$observed)
pdf(out, width = 12, height = height)
ggplot(data = TABLE,aes(x=FDR, 
           y = reorder(name,-1*extra_group_order_new),
           #y=reorder(name,-1*FDR), 
           colour=observed, 
           size=enrichment)) +
  geom_point() +
  scale_fill_gradient(limits = c(0,max), breaks=seq(0,max)) +
  scale_size(range = c(5, 15)) +
  xlim(0, 0.05) + 
  labs(x="FDR", y="", colour="DE genes", size="Enrichment") +
  theme_bw() +
  theme(axis.text.y = element_text(hjust=0))
dev.off()

