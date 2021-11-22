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


#creation of exp matrix
"
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

library(stringi)
for (filename in list.files(path = './raw/')) {
  sampname <- stri_extract(regex = 'GSM[A-Za-z0-9]*', str = filename)
  samptable <- read.table(file = paste0('./raw/',filename), header = T)
  exp[,sampname] <- samptable$NumReads
}

exp.export <- exp
exp.export <- cbind(rownames(exp.export), exp.export)
exp.export <- rbind(colnames(exp.export), exp.export)
colnames(exp.export)[1] <- 'txid'
write.table(x = exp.export, file = 'VitD.ExpMatrixCounts.tsv', sep = '\t', row.names = F,col.names = F)
"

exp <- read.table(file = 'VitD.ExpMatrixCounts.tsv', header = T, row.names = 1)
exp
typeof(exp[,1])

todelete <- as.logical(rep('FALSE', 213995))
all(exp[1,] == 0)

#remove genes that have count = 0 in all of the samples
for (i in c(1:213995)) {
  if (all(exp[i,] == 0)) {
    todelete[i] <- T
  }
}



filtered.exp <- exp[-which(todelete), ]
is.character(exp[,1])
as.numeric(exp[,1])

typeof(exp[,1])
mean(exp[,1])
colMeans(filtered.exp[,1])

library(stringi)
sampnames <- gse$geo_accession
vitlevel <- as.numeric(stri_extract(str = gse@phenoData@data[["characteristics_ch1.1"]],regex = '[:digit:]+'))
sampleinfo <- data.frame(samplenames = sampnames, vitlevel = vitlevel, oldrows = c(1:98))
sampleinfo <- sampleinfo[order(sampleinfo$vitlevel),]
rownames(sampleinfo) <- NULL
sampleinfo$label <- c(rep('low', 41), rep('high', 57))
hist(sampleinfo$vitlevel[1:41])
hist(sampleinfo$vitlevel[42:98])
sampleinfo <- sampleinfo[order(sampleinfo$oldrows),]
sampleinfo
y <- as.factor(sampleinfo$label)


hist(vitlevel)
hist(vitlevel[seq(1,98,2)])
hist(vitlevel[seq(2,98,2)])

mean(vitlevel[seq(1,98,2)])
median(vitlevel[seq(1,98,2)])
mean(vitlevel[seq(2,98,2)])
median(vitlevel[seq(2,98,2)])

boxplot(vitlevel)
median(vitlevel)
mean(vitlevel)



ynum <- y
levels(ynum) <- c(1,0)
ynum <- as.numeric(ynum)
pca <- prcomp(t(filtered.exp))
plot(pca$x[,c(1:2)],col=4-ynum)

library(DESeq2)
exp.matrix <- filtered.exp
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
signaturegenes



library(ensembldb)
txdb <- makeTxDbFromGFF(file = 'gencode.v38.chr_patch_hapl_scaff.annotation.gff3')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
genename.txname <- tx2gene[which(tx2gene$TXNAME %in% rownames(signaturegenes)),]


signaturegenes <- data.frame(signaturegenes)
signaturegenes <- cbind(rownames(signaturegenes), signaturegenes)
colnames(signaturegenes)[1] <- 'txname'
rownames(signaturegenes) <- NULL
signaturegenes

signaturegenes <- merge.data.frame(x = signaturegenes,y = tx2gene,
                                   by.x = "txname", by.y = "TXNAME", all.x = T)
signaturegenes

signaturegenes <- signaturegenes[order(signaturegenes$padj),]
signaturegenes
rownames(signaturegenes) <- NULL
signaturegenes

#signaturegenes <- na.omit(signaturegenes)

genesig.export <- signaturegenes
genesig.export <- rbind(colnames(genesig.export), genesig.export)
rownames(genesig.export) <- NULL
write.table(x = genesig.export,file = 'vitdgenes.lowhigh.deseq2.tsv', col.names = F, row.names = F)


library(randomForest)

rf2 <- randomForest(t(filtered.exp), y = y, ntree=1000)
rf2 #randomforest 
summary(rf2)
plot(rf2)
#randomforest does not appear to give reliable results, especially when it comes to 
#predicting the low-scored vitamin d patients
#one reason might be that in the split there are more vitamin d high patients
#as compared to vitamin d low patients
