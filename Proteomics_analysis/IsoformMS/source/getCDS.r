list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE) 


in.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/IsoformMS/processed/test/SNX2.gff3"
out.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/IsoformMS/processed/test/SNX2.CDS.gff3"
coord="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/GRCh38_p12/nanopore.SQANTI3.proteome.lookup.txt"

in.path=args[1]
out.path=args[2]
coord=args[3]

GFF=read.table(in.path, sep="\t")

COORD=read.table(coord)

transcripts=substring(do.call(rbind,strsplit(GFF$V9[GFF$V3=="transcript"],split=";"))[,1],4)
strand=unique(GFF$V7)


COORD=COORD[COORD$V1 %in% transcripts,]
COORD$V2=COORD$V2-1

for (t in 1:nrow(COORD)){
  transcript=COORD$V1[t]
  start=COORD$V2[t]
  end=COORD$V3[t]  

  EXONS=GFF[GFF$V3=="exon" & grepl(transcript,GFF$V9),]
  EXONS$len=EXONS$V5-EXONS$V4+1
  
  if (strand=="+"){
    EXONS=EXONS[order(EXONS$V4,decreasing = F),]
  }else{
    EXONS=EXONS[order(EXONS$V5,decreasing = T),]
  }
  
  EXONS$cumsum_1=cumsum(c(0,EXONS$len[-nrow(EXONS)]))
  EXONS$cumsum_2=cumsum(EXONS$len)
  
  select=EXONS$cumsum_2>=start & EXONS$cumsum_1<=end
  
  CDS=EXONS[select,]
  CDS$V3="CDS"
  
  if (strand=="+"){
    CDS$V4[which.min(CDS$V4)]=CDS$V4[which.min(CDS$V4)] + (start - CDS$cumsum_1[which.min(CDS$V4)]) 
    
    CDS$V5[which.max(CDS$V5)]=CDS$V5[which.max(CDS$V5)]  - (CDS$cumsum_2[which.max(CDS$V5)]-end ) 
    
    
  }else{ 
    CDS$V5[which.max(CDS$V5)]=CDS$V5[which.max(CDS$V5)] - (start-CDS$cumsum_1[which.max(CDS$V5)]) 
    CDS$V4[which.min(CDS$V4)]=CDS$V4[which.min(CDS$V4)] + (CDS$cumsum_2[which.min(CDS$V4)]-end ) 
  }
  
  CDS$len2=CDS$V5-CDS$V4+1
  
  CDS$cumsum_3=cumsum(CDS$len2)
  
  CDS$modolo=CDS$cumsum_3 %% 3
  
  CDS$V8=c(0,(3-CDS$modolo[-nrow(CDS)])%%3)
  
  if(CDS$modolo[nrow(CDS)]!=0){
    print("ERROR!!!!!!!!!!")
  }
  
  GFF=rbind(GFF,CDS[,1:9])
  
}



write.table(GFF, col.names = F, row.names = F, sep = "\t", quote=F, file=out.path) 

