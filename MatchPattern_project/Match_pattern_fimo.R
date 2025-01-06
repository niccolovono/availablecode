library(memes)
library(universalmotif)
library(IRanges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19 

test_seq <- get_sequence(genome = hg19, regions = "chr4:100206129-100242558", )

res <- runFimo(sequences = test_seq, motifs = "~/Downloads/MA0158.2.meme", parse_genomic_coord = T)
resdf <- as.data.frame(res)
