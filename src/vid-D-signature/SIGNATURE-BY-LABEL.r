'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")'

'if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ensembldb")'

getwd()
setwd("~/Desktop/dmlab/")
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)
library(GEOquery)
#getGEO(GEO = "GSE157982", destdir = './')
gse <- getGEO(filename = './GSE157982_series_matrix.txt.gz')

exp <- data.frame(matrix(ncol = 98, nrow = 213995))
#library(stringi)
#sampnames <- stri_extract(regex = 'GSM[A-Za-z0-9]*', str = list.files(path = './GSE157984_RAW/'))
sampnames <- gse$geo_accession
sampnames

colnames(exp) <- sampnames

tab <- read.table(file = './GSM4783540_MD11535_B.sf.txt', header = T)
tab
rownames(exp) <- tab$Name


#http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat

for (filename in list.files(path = './raw/')) {
  sampname <- stri_extract(regex = 'GSM[A-Za-z0-9]*', str = filename)
  samptable <- read.table(file = paste0('./raw/',filename), header = T)
  exp[,sampname] <- samptable$NumReads
}

exp.export <- exp
exp.export <- cbind(rownames(exp.export), exp.export)
exp.export <- rbind(colnames(exp.export), exp.export)
colnames(exp.export)[1] <- 'txid'
write.table(x = exp.export, file = 'VitD.ExpMatrixCounts.tsv')

y <- as.factor(gse@phenoData@data[["characteristics_ch1"]])
head(y)
levels(y) <- c('post','pre')
head(y)

ynum <- y
levels(ynum) <- c(1,0)
ynum <- as.numeric(ynum)
pca <- prcomp(t(exp))
plot(pca$x[,c(1:2)],col=4-ynum)

library(DESeq2)
exp.matrix <- exp
rownames(exp) <- NULL
colnames(exp) <- NULL


sample.info <- data.frame(condition = y)
rownames(sample.info) <- sampnames
all(rownames(sample.info) == colnames(exp.matrix))
#?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData = round(exp.matrix), colData = sample.info, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
rownames(res)[1]
exp[1,]
length(which(res$padj < 0.05))

signaturegenes <- res[which(res$padj < 0.05),]
signaturegenes <- signaturegenes[order(signaturegenes$padj),]



library(ensembldb)
txdb <- makeTxDbFromGFF(file = 'gencode.v38.chr_patch_hapl_scaff.annotation.gff3')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
genename.txname <- tx2gene[which(tx2gene$TXNAME %in% rownames(signaturegenes)),]


signaturegenes <- data.frame(signaturegenes)
signaturegenes <- cbind(rownames(signaturegenes), signaturegenes)
colnames(signaturegenes)[1] <- 'txname'
rownames(signaturegenes) <- NULL


signaturegenes <- merge.data.frame(x = signaturegenes,y = tx2gene,
                                   by.x = "txname", by.y = "TXNAME", all.x = T)
signaturegenes

signaturegenes <- signaturegenes[order(signaturegenes$padj),]
signaturegenes
rownames(signaturegenes) <- NULL
signaturegenes


notfound <- which(is.na(signaturegenes$GENEID))
signaturegenes[notfound,"txname"]
signaturegenes[notfound,"GENEID"] <- 'ENSG00000196586.16'


genesig.export <- signaturegenes
genesig.export <- rbind(colnames(genesig.export), genesig.export)
colnames(genesig.export) <- NULL

write.table(x = genesig.export,file = 'vitdgenes.deseq2.tsv', col.names = F, row.names = F)

