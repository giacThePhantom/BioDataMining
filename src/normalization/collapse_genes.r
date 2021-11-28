#!/bin/Rscript

library(tidyverse)

read_exp_data <-function(filename){
  res <- read.csv(filename)
  return(res)
}

write_collapsed_genes <- function(filename, df){
  write.csv(df, filename)
}

get_gene_list <- function(df){
  res <- unique(df$Symbols)
  return(res)
}

collapse_gene <- function(df, gene_list){
  res <- df
  for(i in gene_list){
    gene_row <- res[res$Symbols == i,]
    if(nrow(gene_row) > 1){
      print("found")
      new_row <- mapply(mean, gene_row[,-c(1,2,3)])  ## To be changed to statistic used
      new_row <- as.list(c(gene_row[1,c(1,2,3)], new_row))
      res <- setdiff(res, gene_row)
      res <- rbind(res, new_row)
    }
  }
  return(res)
}








args <- commandArgs(trailingOnly = TRUE)

exp_data <- read_exp_data(args[1])

genes <- get_gene_list(exp_data)

collapsed_data <- collapse_gene(exp_data, genes)

write_collapsed_genes(args[2], collapsed_data)
