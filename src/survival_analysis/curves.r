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
		df[,gene_cutoff[i,"genes"]][df[,gene_cutoff[i,"genes"]] > gene_cutoff[i, "cutpoint"]] <- "high"
		df[,gene_cutoff[i,"genes"]][df[,gene_cutoff[i,"genes"]] <= gene_cutoff[i, "cutpoint"]] <- "low"
	}
	return(df)
}

create_plot <- function(formula, data){
	fit <- surv_fit(formula, data = df)
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

save_plot <- function(plot, name){
	name <- paste0(name, ".png")
	png(name)
	print(plot, newpage = FALSE)
	dev.off()


}


args <- commandArgs(trailingOnly = TRUE)

gene_cutoff <- get_gene_cutoff(args[1])

temp <- data.frame(cutpoint = c(0), genes = c("sex"))

#
df <- get_surv_data(args[2], c(gene_cutoff$genes))

df <- split_high_low(df, gene_cutoff)


for(i in temp$genes){
	form <- as.formula(paste("Surv(time, status) ~", i))
	plot <- create_plot(form, df)
	save_plot(plot, paste0(args[3], "/", i))
}
#
#
#name <- "temp"
