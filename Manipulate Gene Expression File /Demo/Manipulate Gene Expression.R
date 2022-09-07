library(tidyverse)
library(DESeq2)
library(airway)

#counts_data describe countable quantities of the genes' expression.
counts_data <- read.csv('~/RScript/R-data/counts_data.csv')

#sample_info represent the design of the study, labeling samples as treated or untreated
colData <- read.csv('~/RScript/R-data/sample_info.csv')

#control if the matches among samples in sample_info and count_data match
all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data)==rownames(colData))

#construct the data-set for further analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data, #array reading data
                       colData = colData,  #characteristics of samples
                       design = ~ dexamethasone) #design of the study

#pre-filtering
keep <- rowSums(counts(dds)) >= 10 #filtering all the samples with reading values > than 10
dds <- dds[keep,] #eliminate all the samples with reading values < than 10

#remember: collapse technical replicates (not biological replicates)
 
#set factor level (treated and untreated)
dds$betamethasone <- relevel(dds$dexamethasone, ref="untreated") #construct the levels of analysis: we want to compare treated samples with untreated samples

#perform DESeq analysis of our filtered and pre-processed data:
dds <- DESeq(dds)

#save the results in a new dataframe - view(res) to see what it looks like
res <- results(dds)

#adjust p-values for significance - just change "alpha" variable.
#tidy = TRUE export the variable as a dataframe
summary(res)
res005 <- results(dds, alpha = 0.05)#, tidy = TRUE) #export a dataframe
summary(res005)


#contrast function
#results(dds, contrast = c("dexamethasone", "leveltobecompared", "untreated"))

#plot a volcano plot with informations about upregulated or downregulated genes retrieved fro mthe analysis
plotMA(res005)


#collect and convert official gene symbols from ensebl id using ensebl database
library(biomaRt)
library(annotables)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

## export the results of the analysis as csv and create a new csv file with the ensebl_ids only

#load and read ensembl ids file
ens_id <- read.csv("ens_id.csv")

#call the ensembl database
listEnsembl()
ensembl <- useEnsembl(biomart = "genes")
datasets <- listDatasets(ensembl)
ensembl.con <- useMart("ensembl", dataset = 'hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
ens_id_conv <- getBM(attributes= c("ensembl_gene_id", "external_gene_name"),
      filters = "ensembl_gene_id",
      values = res005$row, #which values we need to 
      mart = ensembl.con)

#merge in one dataframe the values from the analysis and the resulting gene symbol retrieved from the ensembl database
data_fin <- merge(res005, ens_id_conv, by.x ="row", by.y="ensembl_gene_id")

#clean the data removing NA and duplicated data (TECHNICAL not biological replicates)
data_fin <- na.omit(data_fin)
data_fin <- data_fin[!duplicated(data_fin),]

#divide in two dataframe the upregulated genes and downregulated genes, it could be useful
upregulated <- data_fin %>%
  filter(log2FoldChange>1)
downregulated <- data_fin %>%
  filter(log2FoldChange<1)

library(dplyr)

#append new column based on other column value
data_fin<- data_fin %>%
  mutate(U_D = case_when(log2FoldChange>1 ~ "upregulated",
                         log2FoldChange<1 ~ "downregulated")
  )

#plot the count of downregulated vs. upregulated genes
ggplot(data_fin, aes(x=U_D)) + 
  geom_bar()

#prepare the data for visualizaton setting a specific foldchange treshold
up_filtered <- upregulated %>%
  filter(log2FoldChange>2)
ggplot(up_filtered, aes(x=log2FoldChange, y=external_gene_name)) + 
  geom_point()
#same can be done for downregulated genes

#append progressive id column to the dataframe
up_filtered <- up_filtered %>%
  mutate(id = c(1:130))

#roundplot
label_data <- up_filtered

number_of_bar <- nrow(label_data)
angle <- 90-360 * (label_data$id-0.5)/number_of_bar
label_data$hjust <- ifelse(angle < -90,1,0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

p <- ggplot(up_filtered, aes(x=as.factor(id), y=log2FoldChange+30)) +
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  ylim(-100,120)+
  theme_minimal() + 
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,5), "cm")
  ) +
  coord_polar(start=0) +
  geom_text(data=label_data, aes(x=id, y=log2FoldChange+60, hjust=hjust), label=up_filtered$external_gene_name, color="black", fontface="bold", alpha=0.6, size=2.5, angle=label_data$angle, inherit.aes = FALSE)




