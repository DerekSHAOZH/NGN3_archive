library(Gviz)
library(GenomicRanges)
library(GenomicFeatures)
# #Load data : class = GRanges)
# data(cpgIslands)
# 
# #Plot annotation track ####
# #Annotation track, title ="CpG"
# atrack <- AnnotationTrack(cpgIslands, name = "CpG")
# plotTracks(atrack)
# 
# #Add genome axis track ####
# ## genomic coordinates
# gtrack <- GenomeAxisTrack()
# plotTracks(list(gtrack, atrack))
# 
# #Add chromosome ideogram ####
# #genome : "hg19" 
# gen<-genome(cpgIslands)
# #Chromosme name : "chr7"
# chr <- as.character(unique(seqnames(cpgIslands)))
# #Ideogram track
# itrack <- IdeogramTrack(genome = gen, chromosome = chr)
# plotTracks(list(itrack, gtrack, atrack))
# 
# #Add gene model ####
# #Load data
# data(geneModels)
# head(geneModels)
# 
# #Plot
# grtrack <- GeneRegionTrack(geneModels, genome = gen,
#                            chromosome = chr, name = "Gene Model")
# plotTracks(list(itrack, gtrack, atrack, grtrack))
# 
# #Building GeneRegionTrack objects from gtf files ####
# fn <- '/project/Neurodifferentiation_System/GeneralResources/gencode.v32.primary_assembly.annotation.gtf'
# options(ucscChromosomeNames=FALSE)
# txdb <- makeTxDbFromGFF(fn)
# gm <- GeneRegionTrack(txdb, chromosome='chr12', name = 'test',
#                       transcriptAnnotation = 'transcript')
# plotTracks(gm, from = 25700000, to = 26850000)
# gm@range@elementMetadata$feature
# 
# 
# 
# #Plotting parameters ####
# #Annotation of transcript
# #Change panel and title background color
# grtrack <- GeneRegionTrack(geneModels, genome = gen,
#                            chromosome = chr, name = "Gene Model", 
#                            # transcriptAnnotation = "symbol",
#                            transcriptAnnotation = 'transcript',
#                            background.panel = "#FFFEDB",
#                            background.title = "darkblue")
# plotTracks(list(itrack, gtrack, atrack, grtrack))
# 
# #Zoom the plot ####
# #Use from and to arguments to zooms
# plotTracks(list(itrack, gtrack, atrack, grtrack),
#            from = 26700000, to = 26750000)
# # Use extend.left and extend.right to zoom
# #those arguments are relative to the currently displayed ranges, 
# #and can be used to quickly extend the view on one or both ends of the plot.
# plotTracks(list(itrack, gtrack, atrack, grtrack),
#            extend.left = 0.5, extend.right = 0.5)
# # to drop the bounding borders of the exons and 
# # to have a nice plot
# plotTracks(list(itrack, gtrack, atrack, grtrack),
#            extend.left = 0.5, extend.right = 1000000, col = NULL)
# 
# #Building DataTrack objects from files ####
# bamFile <- system.file("extdata/test.bam", package = "Gviz")
# dTrack4 <- DataTrack(range = bamFile, genome = "hg19",
#                      type = "l", name = "Coverage", window = -1, chromosome = "chr1")
# plotTracks(dTrack4, from = 189990000, to = 190000000)
# 
# #Building AnnotationTrack objects from bam files ####
# #Annotation track
# aTrack2 <- AnnotationTrack(range = bamFile, genome = "hg19",
#                            name = "Reads", chromosome = "chr1")
# plotTracks(aTrack2, from = 189995000, to = 190000000)
# 
# plotTracks(list(dTrack4, aTrack2), from = 189990000,
#            to = 190000000)
# 
# #AlignmentsTrack ####
# afrom=2960000
# ato=3160000
# #bam file
# alTrack <- AlignmentsTrack(system.file(package = "Gviz", "extdata", "gapped.bam"), isPaired = TRUE)
# bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr12",
#                               start = afrom, end = ato, filter = list(with_ox_refseq_mrna = TRUE),
#                               stacking = "dense")
# 
# plotTracks(c(gm, alTrack), from = afrom, to = ato, chromosome = "chr12")
# plotTracks(c(gm, alTrack), from = afrom + 12700,
#            to = afrom + 15200, chromosome = "chr12")
# plotTracks(c(gm, alTrack), from = afrom + 12700, 
#            to = afrom + 15200, chromosome = "chr12",
#            type = c("coverage", "sashimi", 'pileup'))
# introns <- GRanges("chr12", IRanges(start = c(2973662, 2973919),
#                                     end = c(2973848, 2974520)))
# plotTracks(c(gm, alTrack), from = afrom + 12700, to = afrom + 15200, 
#            chromosome = "chr12", type = c("coverage", "sashimi", 'pileup'), 
#            sashimiFilter = introns)
# plotTracks(c(gm, alTrack), from = afrom + 12700, to = afrom + 15200, 
#            chromosome = "chr12", type = c("coverage", "sashimi", 'pileup'), 
#            sashimiFilter = introns, sashimiFilterTolerance = 5L,
#            col.mates = "purple", col.gap = "grey")
# 
# #Track highlighting and overlays ####
# #highlight
# 
# ht <- HighlightTrack(trackList = list(atrack, grtrack,dTrack4),
#                      start = c(26705000, 26720000), width = 7000, chromosome = 'chr7')
# plotTracks(list(itrack, gtrack, ht), chromosome = 'chr7', 
#            from = 26700000, to = 26750000)
# 





#Application: MMP23B day5_1 ####
plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Plots/'
gtf_fn = '/project/Neurodifferentiation_System/GeneralResources/MMP23B.gff3'
options(ucscChromosomeNames=FALSE)
txdb <- makeTxDbFromGFF(gtf_fn)
grtrack <- GeneRegionTrack(txdb, chromosome='chr1', name = 'MMP23B',
                      transcriptAnnotation = 'transcript',
                      min.height =5)


                      
plotTracks(grtrack, from = 1631000, to = 1636000, stackHeight = 0.3)
gtrack <- GenomeAxisTrack()
#something went wrong when try to access hg38 genome assembly
# itrack <- IdeogramTrack(genome = 'hg19', chromosome = 'chr1')
# plotTracks(list(itrack, gtrack, grtrack), from = 1631000, to = 1636000)
plotTracks(list(gtrack, grtrack), from = 1631000, to = 1636000, sizes=c(5,1))

rna_fn <- '/project/owlmayerTemporary/Derek/RNA_alignGenome/alignment/day5_1/Aligned.sortedByCoord.out.bam'
ont_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/IGV/genome/day5_1.nonred.genome.sorted.bam'

intron <- GRanges("chr1", IRanges(start = 1632374,
                                  end = 1632783))

rna_align <- AlignmentsTrack(rna_fn, isPaired = TRUE, chromosome = 'chr1', name = "RNA-seq (day 5)",
                             sashimiFilter = intron, sashimiFilterTolerance = 3L, 
                             lwd.sashimiMax = 1,
                             col.gap = "lightblue1", col.reads = 'gray20'
                             # stackHeight = 0.1#,
                             # type = c('pileup', 'coverage',"sashimi")
                             )
ont_align <- AlignmentsTrack(ont_fn, isPaired = FALSE, chromosome = 'chr1', name = 'ONT-seq (day 5)',
                             sashimiFilter = intron, sashimiFilterTolerance = 3L,
                             # stackHeight = 0.1, 
                             col.gap = "lightblue1", col.reads = 'gray20'#,
                             # type = c('pileup', 'coverage',"sashimi" )
                             )


ht <- HighlightTrack(trackList = list(grtrack, ont_align, rna_align),
                     start = 1632374, width = 409, chromosome = 'chr1',
                     inBackground = T)

pdf(paste0(plot_dir_out,'gviz_MMP23B.pdf'), 6, 5.5)
plotTracks(list(gtrack, ht), from = 1631300, to = 1635000, title.width = 2,
           chromosome = "chr1", sizes = c(0.7,1,3,3), min.height = 1.5,
           # sashimiFilter = intron, #sashimiFilterTolerance = 5L,
           coverageHeight = 0.15, minCoverageHeight = 0,
           sashimiHeight=0.13, minSashimiHeight = 0,
           col.sashimi = 'red', type = c('pileup', 'coverage',"sashimi")
           )
dev.off()



#Application: DHX15 day5_3 ####
plot_dir_out <- '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Plots/'
gtf_fn = '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/SQANTI3_filter/nanopore.filtered.gtf'
options(ucscChromosomeNames=FALSE)
txdb <- makeTxDbFromGFF(gtf_fn)
grtrack <- GeneRegionTrack(txdb, chromosome='chr4', name = 'DHX15',
                           transcriptAnnotation = 'transcript',
                           min.height =5)



plotTracks(grtrack, from = 24525475, to = 24586554, stackHeight = 0.3)
gtrack <- GenomeAxisTrack()
#something went wrong when try to access hg38 genome assembly
# itrack <- IdeogramTrack(genome = 'hg19', chromosome = 'chr1')
# plotTracks(list(itrack, gtrack, grtrack), from = 1631000, to = 1636000)
# plotTracks(list(gtrack, grtrack), from = 1631000, to = 1636000, sizes=c(5,1))

# rna_fn <- '/project/owlmayerTemporary/Derek/RNA_alignGenome/alignment/day5_1/Aligned.sortedByCoord.out.bam'
ont_fn <- '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/IGV/ONT_genome_alignment/day5_3.nonred.genome.sorted.bam'

# intron <- GRanges("chr1", IRanges(start = 1632374,
                                  # end = 1632783))

# rna_align <- AlignmentsTrack(rna_fn, isPaired = TRUE, chromosome = 'chr1', name = "RNA-seq (day 5)",
#                              sashimiFilter = intron, sashimiFilterTolerance = 3L, 
#                              lwd.sashimiMax = 1,
#                              col.gap = "lightblue1", col.reads = 'gray20'
#                              # stackHeight = 0.1#,
#                              # type = c('pileup', 'coverage',"sashimi")
# )
ont_align <- AlignmentsTrack(ont_fn, isPaired = FALSE, chromosome = 'chr1', name = 'ONT-seq (day 5)',
                             # sashimiFilter = intron, sashimiFilterTolerance = 3L,
                             # stackHeight = 0.1,
                             stacking = "squish",
                             col.gap = "lightblue1", col.reads = 'gray20'#,
                             # type = c('pileup', 'coverage',"sashimi" )
)


ht <- HighlightTrack(trackList = list(grtrack, ont_align),
                     start = 24537048, width = 3251, chromosome = 'chr4',
                     inBackground = T)

pdf(paste0(plot_dir_out,'gviz_DHX15.pdf'), 6, 5)
plotTracks(list(gtrack, ht), from = 24525475, to = 24586554, title.width = 2,
           chromosome = "chr4", sizes = c(1,1,3), min.height = 1.5,
           # sashimiFilter = intron, #sashimiFilterTolerance = 5L,
           coverageHeight = 0.15, minCoverageHeight = 0,
           # sashimiHeight=0.13, minSashimiHeight = 0,
           # col.sashimi = 'red', 
           # type = c('pileup', 'coverage',"sashimi")
           type = c('pileup', 'coverage')
)
dev.off()

