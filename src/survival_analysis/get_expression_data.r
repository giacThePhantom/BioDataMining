#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

get_exp_levels <- function(filename){
	res <- read.csv(filename, header = FALSE)
	res <- as.data.frame(t(res))
	res <- type.convert(res, as.is = TRUE)
	res[4,1] <- "sample"
	colnames(res) <- unlist(res[4,])
	res <- res[-c(1,2,3,4),]
	return(res)
}

merge_exp_clin_data <- function(exp_data, filename){
	res <- read.csv(filename)
	res <- merge(exp_data, res, by.x="sample", by.y="sample", all = TRUE)
	res <- type.convert(res, as.is = TRUE)
	return(res)

}

exp_data <- get_exp_levels(args[1])
	
write.csv(merge_exp_clin_data( exp_data, args[2]), args[3])
