library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db) # For mapping gene symbols to Entrez IDs

# Load the TxDb object
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Extract gene information
genes <- genes(txdb)

# View the gene data
head(genes)

# Get Entrez IDs for gene symbols
symbols <- c("ADH1B") # Example symbols
entrez_ids <- mapIds(org.Hs.eg.db, keys = symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Subset the GRanges object by Entrez IDs
gene_locations <- genes[names(genes) %in% entrez_ids]

# Add gene symbols as metadata
mcols(gene_locations)$symbol <- symbols[match(names(gene_locations), entrez_ids)]
gene_locations

gene_locations@ranges
motif_file <- "~/Downloads/MA0158.2.meme"

#----

adjust_start <- 20000  # Number to subtract from start position

# Prepare sequences and run FIMO
results_list <- lapply(seq_along(gene_locations), function(i) {
  # Extract gene information
  gene <- gene_locations[i]
  chr <- as.character(seqnames(gene))
  start <- start(gene) - adjust_start
  end <- end(gene)
  
  # Ensure start is not negative
  start <- ifelse(start < 0, 0, start)
  
  # Construct the regions argument
  regions <- sprintf("%s:%d-%d", chr, start, end)
  
  # Get sequence
  test_seq <- get_sequence(genome = hg19, regions = regions)
  
  # Run FIMO
  res <- runFimo(sequences = test_seq, motifs = motif_file, parse_genomic_coord = TRUE)
  
  # Convert results to a data frame
  as.data.frame(res)
})

# Combine all results into a single data frame
results_df <- bind_rows(results_list)
