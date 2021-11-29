library(GEOquery)
library(affy)
library(frma)
library("hgu133afrmavecs")

setwd("//percorso//file//CEL")

f <- list.files(pattern = ".CEL")

i = 0
j = 0

for(i = 1, i <= j, j++ in 1:number of file){ # set number of file in the folder
  affy.data <- ReadAffy(verbose = FALSE, filenames = f[i])
  frma_data <- frma(affy.data)
  exprs <- assayData(frma_data)$exprs
  }

write.table(exprs, file = "frma.txt", quote = FALSE, sep = "\t")

//Annotation

probes=row.names(frma)
ls("package:hgu133a.db")
Symbols = unlist(mget(probes, hgu133aSYMBOL, ifnotfound=NA))
Entrez_IDs = unlist(mget(probes, hgu133aENTREZID, ifnotfound = NA))
frma = cbind(probes, Symbols, Entrez_IDs, frma)
write.table(frma, file = "annotation_frma.txt", sep ="\t", row.names = FALSE, col.names = TRUE)
