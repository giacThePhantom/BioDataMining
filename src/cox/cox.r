#!/bin/Rscript

library(survival)
library(survminer)
library(tidyverse)

############################ TO DO ###############################
#																															   #
#  1. Use both datasets																					 #
#  2. Collapse gene (see if Gaia can do it better or use median) #
#  3. Filter data																								 #
#  4. Transform as numeric																			 #
#  5. Create the cox function with data													 #
#																																 #
##################################################################



get_clin_data <- function(filename, dataset, col_filter){
	samples_clinical_data <- read.csv(filename, header = TRUE)
	print(colnames(samples_clinical_data))
	samples_clinical_data <- samples_clinical_data[samples_clinical_data$dataset %in% dataset,]
	return(samples_clinical_data[,col_filter])
}




get_exp_data <- function(filename){
	rma_data <- read.csv(filename, header = FALSE)
	res <- rma_data[complete.cases(rma_data[2]), ]
	res <- as.data.frame(t(res))
	res[2,1] <- "sample"
	colnames(res) <- res[2,]
	res <- res[-c(1,2,3),]

	print(tibble(res))
	res <- transform(res, DDR1 = as.numeric(DDR1))
	#Collapse genes
	return(res)
}

merge_exp_clin <- function(exp_data, clin_data){
	res <- merge(exp_data, clin_data, by.x="sample", by.y="geo_accession", all = TRUE)
	res <- transform(res, overall_surv_months = as.numeric(overall_surv_months))
	res <- transform(res, overall_surv_months = as.numeric(overall_surv))
	return(res)
}


args <- commandArgs(trailingOnly=TRUE)


dataset <- c("GSE17536")
col_filter <- c("sex",
								"stage",
								"overall_surv",
								"dis_free_surv",
								"dis_specific_surv",
								"geo_accession",
								"overall_surv_months",
								"dis_free_surv_months",
								"dis_specific_surv_months"
								)

clin_data <- get_clin_data(args[1], dataset, col_filter)

exp_data <- get_exp_data(args[2])

cox_data <- merge_exp_clin(exp_data, clin_data)

print(tibble(clin_data))
print(tibble(exp_data))
print(tibble(cox_data))

write.csv(cox_data, file="test.csv" )

cox_data <- cox_data[complete.cases(cox_data[, c("overall_surv", "overall_surv_months", "DDR1")]),]
print(tibble(cox_data))
#print(head(samples_clinical_data))

cox <- coxph(Surv(overall_surv_months, overall_surv) ~ DDR1, data = cox_data)

print(cox)


#covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")


#univ_formulas <- sapply(covariates,
#                        function(x) as.formula(paste('Surv(time, status)~', x)))



#univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})



#univ_results <- lapply(univ_models,
#                       function(x){
#                          x <- summary(x)
#                          p.value<-signif(x$wald["pvalue"], digits=2)
#                          wald.test<-signif(x$wald["test"], digits=2)
#                          beta<-signif(x$coef[1], digits=2);#coeficient beta
#                          HR <-signif(x$coef[2], digits=2);#exp(beta)
#                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
#                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
#                          HR <- paste0(HR, " (",
#                                       HR.confint.lower, "-", HR.confint.upper, ")")
#                          res<-c(beta, HR, wald.test, p.value)
#                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
#                                        "p.value")
#                          return(res)
#                          #return(exp(cbind(coef(x),confint(x))))
#                         })
#res <- t(as.data.frame(univ_results, check.names = FALSE))
#as.data.frame(res)
