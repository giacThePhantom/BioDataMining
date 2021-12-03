#!/bin/Rscript

library(survival)
library(survminer)
library(survMisc)
library(plyr)
library(tidyverse)
library(readr)
library(sos)
library(scales)
library(shiny)
library(`ggplot2`)
library(plotly)

get_clin_data <- function(filename, dataset, col_filter){
	res <- ldply(filename, read_csv)
	res <- res[res$dataset %in% dataset,]
	res <- res[,col_filter]
	return(res)
}




get_exp_data <- function(filename){
	res <- NULL
	for(i in filename){
		rma_data <- read.csv(i, header = FALSE)
		rma_data <- rma_data[complete.cases(rma_data[3]), ]
		rma_data <- as.data.frame(t(rma_data))
		rma_data <- type.convert(rma_data, as.is = TRUE)
		rma_data[3,1] <- "sample"
		colnames(rma_data) <- as.vector(rma_data[3,])
		rma_data <- rma_data[-c(1,2,3),]
		if(is.null(res)){
			res <- rma_data
		}
		else{
			res <- rbind(res, rma_data)
		}
	}
	#Collapse genes
	return(res)
}

merge_exp_clin <- function(exp_data, clin_data){
	res <- merge(exp_data, clin_data, by.x="sample", by.y="geo_accession", all = TRUE)
	res <- type.convert(res, as.is = TRUE)
	return(res)
}

univariate_cox <- function(cox_data, covariates){
	univ_formulas <- sapply(covariates,
                        	function(x) as.formula(paste('Surv(overall_surv_months, overall_surv)~', x)))
	univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_data)})
	univ_results <- lapply(univ_models,
                       	function(x){
														result.frame <- data.frame(cutp(x))
														colnames(result.frame) <- lapply(colnames(result.frame), function(x) paste0("`", x, "`"))
														cutpoint <- signif( result.frame[[ 1, 1 ]], digits = 4 )
                          	x <- summary(x)
                          	p.value<-signif(x$wald["pvalue"], digits=2)
                          	wald.test<-signif(x$wald["test"], digits=2)
                          	beta<-signif(x$coef[1], digits=2);#coeficient beta
                          	HR <-signif(x$coef[2], digits=2);#exp(beta)
                          	HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          	HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          	HR <- paste0(HR, " (",
                                       	HR.confint.lower, "-", HR.confint.upper, ")")
                          	res<-c(beta, HR, wald.test, p.value, cutpoint)
                          	names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                        	"p.value", "cutpoint")
                          	return(res)
                         	})
	res <- t(as.data.frame(univ_results, check.names = FALSE))
	res <- as.data.frame(res)
	rownames(res) <- rep(1:nrow(res))
	res$genes <- lapply(covariates, function(x) str_remove_all(x, "`"))
	res$genes <- lapply(res$genes, function(x) gsub("-", "\\.", x))
	res$genes <- unlist(res$genes)
	return(res)
}

multivariate_cox <- function(cox_data, multivar){
	formula <- "Surv(overall_surv_months, overall_surv)~"
	for(i in multivar){
		formula <- paste(formula, i, '+', sep='')
	}
	formula <- str_sub(formula, 1, nchar(formula) -1)

	res <- coxph(as.formula(formula), data = cox_data)
	return(res)

}

get_parameters <- function(filename){
	parameters <- read.csv(filename, header = TRUE)
	res <- as.list(parameters)
	res <- lapply(res, function(x){x[nzchar(x)]})
	res <- lapply(res, function(x){na.omit(x)})

	return(res)
}

args <- commandArgs(trailingOnly=TRUE)

parameters <- get_parameters(args[1])
exp_data <- get_exp_data(parameters$exp_data)
clin_data <- get_clin_data(parameters$clin_data, parameters$dataset, parameters$clin_data_col_filter)


cox_data <- merge_exp_clin(exp_data, clin_data)

parameters$genes <- intersect(colnames(cox_data), parameters$genes)
parameters$genes <- lapply(parameters$genes, function(x) gsub("-", ".",  x))
colnames(cox_data) <- lapply(colnames(cox_data), function(x) gsub("-", ".",  x))

cox_values <- univariate_cox(cox_data, parameters$genes)
multi_cox <- multivariate_cox(cox_data, parameters$genes)
print(cox_values)
write.csv(cox_values, args[2])
