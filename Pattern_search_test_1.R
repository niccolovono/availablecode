# Load required libraries
library(GenomicFeatures)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)

# Input: List of gene symbols and patterns
up_genes <- read.csv(file = "~/Downloads/GENES_DOWN.csv", header = T)
colnames(up_genes) <- "genes"

converted_gene_names <- genome_genes %>% 
  dplyr::inner_join(up_genes, by = c("symbol" = "genes")) %>% 
  dplyr::select(symbol, ensgene, entrez, chr, start, end, strand, description)

gene_symbols <- converted_gene_names
patterns <- c("GTAATTAC") # Replace with your patterns

# Create TxDb object for the human genome (GRCh37/hg19)
txdb <- txdbmaker::makeTxDbFromUCSC(genome = "hg19", tablename = "ensGene")

# Get the transcripts for the specified genes
gene_tx <- genes(txdb)
gene_tx <- gene_tx[gene_tx$gene_id %in% converted_gene_names$ensgene]

# Define upstream range (20 kb)
upstream_range <- 20000

# Extract sequences and find matches
results <- tibble(Gene = character(),
                  Chromosome = character(),
                  Start = integer(),
                  End = integer(),
                  Strand = character(),
                  Pattern = character(),
                  MatchStart = integer(),
                  MatchEnd = integer())

for (gene in converted_gene_names$ensgene) {
  # Filter gene-specific information
  gene_info <- gene_tx[gene_tx$gene_id == gene]
  
  for (i in seq_along(gene_info)) {
    chr <- as.character(seqnames(gene_info[i]))
    start <- start(gene_info[i])
    end <- end(gene_info[i])
    strand <- as.character(strand(gene_info[i]))
    
    # Calculate upstream coordinates based on strand
    if (strand == "+") {
      upstream_start <- max(1, start - upstream_range)
      upstream_end <- start - 1
    } else {
      upstream_start <- end + 1
      upstream_end <- end + upstream_range
    }
    
    # Extract sequence
    upstream_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, chr, upstream_start, upstream_end)
    
    # Match patterns
    for (pattern in patterns) {
      matches <- matchPattern(pattern, upstream_seq)
      if (length(matches) > 0) {
        match_data <- tibble(
          Gene = gene,
          Chromosome = chr,
          Start = upstream_start,
          End = upstream_end,
          Strand = strand,
          Pattern = pattern,
          MatchStart = start(matches),
          MatchEnd = end(matches)
        )
        results <- bind_rows(results, match_data)
      }
    }
  }
}

1# Save or return results
print(results)
write.csv(results, "gene_pattern_matches.csv", row.names = FALSE)
