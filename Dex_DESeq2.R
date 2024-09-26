library(DESeq2)
library(tidyverse)

#installing airway package from BiocManager
BiocManager::install("airway")
library(airway)

#dowloading Data files from airway. code to getData obtained from Github(bioinformagician)
data(airway)
airway
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[,c(2,3)]
sample_info$dex <- gsub('trt', 'treated', sample_info$dex)
sample_info$dex <- gsub('untrt', 'untreated', sample_info$dex)
names(sample_info) <- c('cellLine', 'dexamethasone')
write.table(sample_info, file = "sample_info.csv", sep = ',', col.names = T, row.names = T, quote = F)

countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ',', col.names = T, row.names = T, quote = F)

#preparing count data................

#read in counts data
counts.data <- read.csv('counts_data.csv')
head(counts.data)

#read sample info to identify samples & genes

colData <- read.csv('sample_info.csv')

#matching row names in colData with colums in counts.data
all(colnames(counts.data) %in% rownames(colData))

#Create a DESeqDataset object..................
dds <- DESeqDataSetFromMatrix(countData = counts.data,
                       colData = colData,
                       design = ~dexamethasone)

dds

#prefiltering: remove rows with low gene counts
#keeping rows that have at least 10 reads total

keep <- rowSums(counts(dds)) >= 10 
#view keep to check genes with false attached
#remove genes with false
dds<- dds[keep, ] 
dds

#set the factor level..tells DESeq2 to compare treated samples with untreated
dds$dexamethasone <- relevel(dds$dexamethasone, ref = 'untreated')
dds$dexamethasone

#Collapse Technical replicates if any
#Running DESeq......................
dds<- DESeq(dds)

#explore results
results <- results(dds)
results 
summary(results)

results0.01<- results(dds, alpha= 0.01)
summary(results0.01)

#visualize with MA plot
plotMA(results)
