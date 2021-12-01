library(GEOquery)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(repr)


my1 <- "GSE17536"
gse1 <- getGEO(my1)
length(gse1)

gse1 <- gse1[[1]]
gse1

pData(gse1)
fData(gse1)
exprs(gse1)

summary(exprs(gse1))

exprs(gse1) <- log2(exprs(gse1))

sampleInfo1 <- pData(gse1)

sampleInfo1 <- rename(sampleInfo1,group = source_name_ch1, patient=characteristics_ch1.1)

corMatrix1 <- cor(exprs(gse1),use="c")
pheatmap(corMatrix1) 

pca1<- prcomp(t(exprs(gse1)))

PCA1 <- cbind(sampleInfo1, pca1$x) %>% 
ggplot(aes(x = PC1, y=PC2, col=group)) + geom_point() 

my2 <- "GSE17537"
gse2 <- getGEO(my2)
length(gse2)

gse2 <- gse2[[1]]
gse2

pData(gse2)
fData(gse2)
exprs(gse2)

summary(exprs(gse2))

exprs(gse2) <- log2(exprs(gse2))

sampleInfo2 <- pData(gse2)

sampleInfo2 <- rename(sampleInfo2,group = source_name_ch1, patient=characteristics_ch1.1)

corMatrix2 <- cor(exprs(gse2),use="c")
pheatmap(corMatrix2) 

pca2<- prcomp(t(exprs(gse2)))

PCA2 <- cbind(sampleInfo2, pca2$x) %>% 
ggplot(aes(x = PC1, y=PC2, col=group)) + geom_point() 





