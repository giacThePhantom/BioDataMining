#!/usr/bin/Rscript

library('biomaRt')

ens_list <- function(filename){
  df <- read.csv(filename, sep = '\t')
  res <- sapply(strsplit(df$GENEID, '\\.'), '[[', 1)
  res <- res[!is.na(res)]
  return(res)
}

mart_init <- function(){
  res <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  return(res)
}

get_genes_symbol <- function(mart, genes){
  filters <- c("ensembl_gene_id")
  attributes <- c("ensembl_gene_id", "description", "hgnc_symbol")
  res <- getBM(filters = filters, attributes = attributes, values = genes, mart = mart)
  return(res)
}


args <- commandArgs(trailingOnly = TRUE)

genes <- ens_list(args[1])

mart <- mart_init()

res <- get_genes_symbol(mart, genes)

print(res)
