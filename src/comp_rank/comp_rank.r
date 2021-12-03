#!/usr/bin/Rscript


get_ranking <- function(filename){
	res <- read.csv(filename)
	return(res)
}


get_signature <- function(filename){
	res <- read.csv(filename)
	res <- res$Symbol
	return(res)
}

args <- commandArgs(trailingOnly = TRUE)

rank <- get_ranking(args[1])

sig <- get_signature(args[2])

pos <- match(sig, rank$Symbol)

avg_rank <- sum(na.omit(pos))/length(pos)

print(avg_rank)
print(nrow(rank)/2)

