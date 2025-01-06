##### ULTIMO DEFINITIVO ######

library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(memes)
library(universalmotif)
library(IRanges)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb)

### I/O
gene_list_path <- "~/Downloads/GENES_DOWN.csv"
output_path_name <- "~/Desktop/"
motif_file <- "~/Downloads/MA0158.2.meme"
adjust_start <- 20000 ## upstream cutoff

## read and convert gene symbols
gene_list <- read.csv(gene_list_path, header = F)
symbols <- gene_list$V1 # Example symbols
entrez_ids <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

#get gene locations
gene_locations <- genes[names(genes) %in% entrez_ids]
mcols(gene_locations)$symbol <- symbols[match(names(gene_locations), entrez_ids)]
gene_locations
gene_locations@ranges

results_list <- lapply(seq_along(gene_locations), function(i) {
  # Extract gene information
  gene <- gene_locations[i]
  chr <- as.character(seqnames(gene))
  start <- start(gene) - adjust_start
  end <- end(gene)
  gene_name <- mcols(gene)$symbol  # Retrieve gene name or ID
  
  # Ensure start is not negative
  start <- ifelse(start < 0, 0, start)
  
  # Construct the regions argument
  regions <- sprintf("%s:%d-%d", chr, start, end)
  
  # Get sequence
  test_seq <- get_sequence(genome = hg19, regions = regions)
  
  # Run FIMO
  res <- runFimo(sequences = test_seq, motifs = motif_file, parse_genomic_coord = TRUE)
  
  # Convert results to a data frame
  if (length(res) > 0) {
    # If matches are detected, add the gene name/ID as a column
    resdf <- as.data.frame(res)
    resdf$gene_name <- gene_name  # Add gene name/ID as a column
  } else {
    # If no matches are detected, return an empty data frame with the gene name column
    resdf <- data.frame(gene_name = gene_name, stringsAsFactors = FALSE)
  }
  
  resdf
})

# Combine all results into a single data frame
results_df <- bind_rows(results_list)

# View the combined results
head(results_df)

write.csv(results_df, file = paste0(output_path_name, "Genes_down_20k.csv"))
