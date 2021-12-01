#!/usr/bin/Rscript

library(votesys)

open_dir <- function(dir){
	res <- list.files(dir)
	res <- lapply(res, function(x) paste(dir, x, sep = '/'))
	return(res)
}


get_raw_data <- function(dir){
	files <- open_dir(dir)
	data <- lapply(files, read.csv)
	data <- lapply(data, function(x) x$Symbols)
	print(files)
	print(head(data[[1]]))
	res <- create_vote(data, xtype = 3, candidate = data[[1]])
	return(res)
}


get_borda_count <- function(vote){
	res <- borda_method(vote, modified = TRUE)
	res <- as.data.frame(res$other_info$count_max)
	res$t <- rownames(res)
	colnames(res) <- c("rank", "Symbol")
	res <- res[order(-res$rank),]
	return(res)
}



args <- commandArgs(trailingOnly = TRUE)

vote <- get_raw_data(args[1])

borda <- get_borda_count(vote)

write.csv(borda, args[2])
