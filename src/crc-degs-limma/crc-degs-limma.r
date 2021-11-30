#!/usr/bin/Rscript

library(tidyverse)
library(limma)


get_clin_data <- function(filename, num){
	res <- read.csv(filename)
	res$stage <- num
	return(res)
}
get_exp_data <- function(filename){
	exp_data <- read.csv(filename)
	exp_data <- type.convert(exp_data, as.is = TRUE)
	colnames(exp_data) <- lapply(colnames(exp_data), function(x) strsplit(x, '_')[[1]][1])
	rownames(exp_data) <- exp_data$Symbols
	return(exp_data)
}


split_exp_data <- function(exp_data, sample){
	sample <- intersect(c('Symbols', sample), colnames(exp_data))
	res <- exp_data[sample]
	return(res)
}


get_samples <-function(low_name, high_name){
	sample_low <- get_clin_data(low_name, 0)
	sample_high <- get_clin_data(high_name, 1)
	res <- rbind(sample_low, sample_high)
	return(res)
}

write_limma <- function(exp, samples, out_name){
	design <- model.matrix(~samples$stage, data = exp)
	fit <- lmFit(exp, design)
	fit <- eBayes(fit)
	fit$Symbols <- rownames(fit)
	fit <- as.data.frame(fit)
	fit <- fit[, c("Symbols", "p.value.samples.stage")]
	colnames(fit) <- c("Symbols", "pvalues")
	fit <- fit[order(fit$pvalues),]
	write.csv(fit, out_name)
}

args <- commandArgs(trailingOnly = TRUE)

samples <- get_samples(args[1], args[2])


exp <- get_exp_data(args[3])

new_samples <- intersect(samples$geo_accession, colnames(exp))
exp <- exp[new_samples]

exp <- exp[, order(names(exp))]

samples <- samples[samples$geo_accession %in% new_samples,]

samples <- samples[order(samples$geo_accession),]

write_limma(exp, samples, args[4])
