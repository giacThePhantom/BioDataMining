setwd('~/Desktop/dmlab/')

#create a single, merged metadata file from the splitted files
"stage.high.metadata <- 'StageHighLow/stage_high.csv'
stage.low.metadata <- 'StageHighLow/stage_low.csv'


metadata.high <- read.csv(file = stage.high.metadata,header = T, sep = ',')
metadata.low <- read.csv(file = stage.low.metadata, header = T,sep = ',')

metadata <- rbind(metadata.high, metadata.low)
metadata

metadata$stage <- c(rep('high', 908), rep('low',2028-908))
metadata$stage[908]
metadata$stage[909]
metadata$stage[2028]

metadata.export <- rbind(colnames(metadata), metadata)
metadata.merged.export <- write.table(x = metadata.export, file = 'metadata.merged.tsv',
                                      row.names = F, col.names = F)

"

#read the complete metadata file (high+low)
input_metadata <- 'metadata.merged.tsv'
metadata <- read.table(file = input_metadata, header = T)

#load the Nth dataset split (eg first bootstrapped dataset)
input_split = 'normalization-final/normalization/stage/annotation/annotation1/annotation_mix1_stage.txt'
split <- read.table(file = input_split, header = T, row.names = 1)

#extract the expression data (colnames=samples, rownames =genenames)
eset <- split[,3:302]

head(eset)
colnames(eset)

#extract sample names from column names
library(stringi)
samplenames <- stri_extract(str = colnames(eset), regex = 'GSM\\d*')
samplenames %in% metadata$geo_accession #check if all sample names are in metadata

which(metadata$geo_accession %in% samplenames) #indices of split samples in metadata excel

split.specific.metadata <- metadata[which(metadata$geo_accession %in% samplenames),]
View(split.specific.metadata)

rownames(split.specific.metadata) <- split.specific.metadata$geo_accession

split.specific.metadata<-split.specific.metadata[samplenames, ] #order metadata
#based on the sequence of samples in expression set

split.specific.metadata$geo_accession == samplenames #check order


#limma workflow
design <- model.matrix(~split.specific.metadata$stage, data = eset)
head(design)

library(limma)
#fit the model
fit <- lmFit(eset, design)
#calculate the t-statistics
fit <- eBayes(fit)
#summarize the results
results <- decideTests(fit)
results
summary(results)


pval <- fit$p.value # probe id with p value
head(pval)
hist(pval)

#add gene names and entrez ID to p value list
filtered.split <- split[,c('Symbols', 'Entrez_IDs')]
head(filtered.split)
split.list <- merge(x = pval, y = filtered.split, by=0, all=T)
colnames(split.list)[c(1,3)] <- c('probeid','pvalue')


#add ensembl gene ID
ensemblid <- read.csv(file = 'mart_export.txt',sep = '\t', header = T)
split.list <- merge(x = split.list, y = ensemblid, by.x = 'Symbols', by.y = 'Gene.name')
head(split.list)
split.list <- split.list[order(split.list$pvalue),] #order list in decreasing p value
rownames(split.list) <- NULL
head(split.list)

#split.list is the final output


"
#additional tests going on..........

#find DEGS using ttests

library('genefilter'')
eset.matrix <- data.matrix(eset)
rownames(eset.matrix) <- NULL
colnames(eset.matrix) <- NULL
eset.matrix <- as.numeric(eset.matrix)
tt40 <- rowttests(x = log(eset.matrix), fac = as.factor(split.specific.metadata$stage))
sum(tt40$p.value < 0.05) #~degs10k
sum(p.adjust(tt40$p.value) < 0.05) #~500degs correcting for the multiple testing

#randomforest variable importance
library(randomForest)
rf <- randomForest(x = t(eset), y = as.factor(split.specific.metadata$stage), ntree = 1000)
rf
plot(rf)
impo <- rf$importance[which(rf$importance != 0),]
impo <- impo[order(impo)]
head(impo)

#lda coefficient estimates
#............."
