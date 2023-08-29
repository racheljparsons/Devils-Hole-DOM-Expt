library(tidyverse)
library(phyloseq)
library(readr)
library(dplyr)
library(ggplot2)


## Removed duplicate samples and Negs. Inverted DH56 (OSU moved samples so 1m was 24m etc) can also remove. Sample DH6 was labeled DH2-12 but that sample was not taken so relabled as DH2-5. Can also remove.
sample.qc <- as.data.frame(read_csv('DH36_meta.csv'))
taxa <- as.data.frame(read_csv('DH36_taxa.csv'))
table.qc <- as.data.frame(read_csv('DH36_table.csv'))
sample.qc$Depth <- as.factor(sample.qc$Depth)

## define the row names for the taxonomy table
row.names(table.qc) <- table.qc$Seq
row.names(taxa) <- taxa$Seq
row.names(sample.qc) <- sample.qc$Seq_ID

## remove the column Seq and ID
table.qc <- table.qc %>% select (-Seq)
taxa <- taxa %>% select (-Seq)
sample.qc <- sample.qc %>% select (-Seq_ID)

## make into matrices
table.qc <- as.matrix(table.qc)
taxa <- as.matrix(taxa)
sample.qc <- as.data.frame(sample.qc)

## transform to phyloseq objects
table.qc <- otu_table(table.qc, taxa_are_rows = TRUE)
taxa <- tax_table(taxa)
sample.qc <- sample_data(sample.qc)

DH36 <- phyloseq(table.qc, taxa, sample.qc)
saveRDS(DH36, file = "DH36.rds")

