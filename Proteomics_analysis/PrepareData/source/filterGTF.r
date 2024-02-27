
filter.gtf.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/GRCh38_p12/nanopore.SQANTI3.gtf"
gtf.path="/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Data/GRCh38_p12/nanopore.new.noSingleExon.gtf"

suppressPackageStartupMessages(library("rtracklayer"))

args <- commandArgs(trailingOnly = TRUE) 


filter.gtf.path=args[1]
gtf.path=args[2]
out=args[3]

                               
filter.gtf=as.data.frame(import(filter.gtf.path, "gtf"))
gtf=import(gtf.path, "gtf")

gtf=gtf[gtf$transcript_id %in% unique(filter.gtf$transcript_id),]


export(gtf, out, "gtf")