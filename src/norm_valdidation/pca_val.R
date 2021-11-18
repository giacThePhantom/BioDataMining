#!/usr/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)

library(tidyverse)

data_files <- list.files(args[1])

print(data_files)

sample_datapoints <- c()

for(f in data_files){
	dataset <- read.csv(paste(args[1], '/', f, sep=''))
	dataset <- dataset[-c(1)]
	print(f)
	temp <- prcomp(dataset)
	print(head(temp$x[1]))
	append(sample_datapoints, temp)
}

print(sample_datapoints)
