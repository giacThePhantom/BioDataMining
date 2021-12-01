#!/bin/Rscript

#clin_data, exp_data, dataset, clin_data_col_filter, genes




clin_data <- c("../../data/splitted_samples/os_data.csv")

exp_data <- c("../../data/collapsed_norm_data/annotation_os_frma.txt")

clin_data_col_filter <- c("stage",
													"overall_surv",
													"dis_free_surv",
													"dis_specific_surv",
													"geo_accession",
													"overall_surv_months",
													"dis_free_surv_months",
													"dis_specific_surv_months")

dataset <- c("GSE17536", "GSE17537")

get_genes <- function(filename, num = -1){
	res <- read.csv(filename)
	res <- res$Symbol
	if(num > 0){
		res <- head(res, n=num)
	}
	return(res)
}

pad_vec <- function(vec, new_length){
	res <- c(vec, rep(NA, new_length-length(vec)))
	return(res)
}

args <- commandArgs(trailingOnly = TRUE)
num = -1
if(length(args) > 2){
	num = as.numeric(args[3])
}

genes <- get_genes(args[1], num)
clin_data <- pad_vec(clin_data, length(genes))
exp_data <- pad_vec(exp_data, length(genes))
clin_data_col_filter <- pad_vec(clin_data_col_filter, length(genes))
clin_data <- pad_vec(clin_data, length(genes))
dataset <- pad_vec(dataset, length(genes))

res <- data.frame(clin_data, exp_data, dataset, clin_data_col_filter, genes)

colnames(res) <- c("clin_data",
									 "exp_data",
									 "dataset",
									 "clin_data_col_filter",
									 "genes")

write.csv(res, args[2])
