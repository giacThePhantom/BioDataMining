#!/usr/bin/Rscript

library(survminer)
library(survival)
library(tidyverse)
library(ggplot2)


get_gene_cutoff <- function(filename){
	res <- read.csv(filename)
	res <- res[res$p.value <= 0.05,]
	res <- res[c("cutpoint", "genes")]
	res$genes = unlist(lapply(res$genes, function(x) gsub("-", ".", x)))
	return(res)
}


get_surv_data <- function(filename, cols){
	res <- read.csv(filename)
	colnames(res) <- unlist(lapply(colnames(res), function(x) gsub("-", ".", x)))
	res <- res[cols]
	return(res)
}

split_high_low <- function(df, gene_cutoff){
	for(i in 1:nrow(gene_cutoff)){
		print(gene_cutoff[i,])
		df[df[gene_cutoff[i,"genes"]] > gene_cutoff[i, "cutpoint"],] <- "high"
		df[df[gene_cutoff[i,"genes"]] <= gene_cutoff[i, "cutpoint"],] <- "low"
	}
	print(df)

}

create_plot <- function(formula, data){
	fit <- survfit(Surv(time, status) ~ sex, data = df)
	res <- ggsurvplot(fit, 
		   	pval = TRUE, 
		   	conf.int = TRUE,
		   	risk.table = TRUE,
		   	risk.table.col = "strata",
		   	surv.median.line = "hv",
		   	ggtheme = theme_bw(),
		   	palette = c("#E7B800", "#2E9FDF"))
	return(res)
}

args <- commandArgs(trailingOnly = TRUE)

gene_cutoff <- get_gene_cutoff(args[1])


print(tibble(gene_cutoff))

temp <- data.frame(cutpoint = c(0), genes = c("sex"))

#
df = read.csv("test.csv")

split_high_low(df, temp)

#
#
#name <- "temp"
#ggsave(device = "pdf", print(plot, newpage = FALSE))
