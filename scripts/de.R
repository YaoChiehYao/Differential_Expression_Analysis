#!/usr/bin/env Rscript
# de.R

library(knitr)
library(readr)
library(tximport)
library(DESeq2)
library(dplyr)
library(tibble)

# Define constants
TESTING <- FALSE # Change to FALSE if using entire Samples set
RESULTS_DIR <- "/home/yao.yao-/BINF6309/module04-YaoChiehYao/results"
AIPTASIA_DIR <- "/work/courses/BINF6309/AiptasiaMiSeq"

# for testing purposes - alternative samples table
testing_samples <- data.frame(Sample = c("Aip02", "Aip02", "Aip02", "Aip02"),
                              Menthol = c("Control", "Control", "Menthol", "Menthol"),
                              Vibrio = c("Control", "Vibrio", "Control", "Vibrio"))
head(testing_samples)

# True script begins
tx2gene <- read.csv(file.path(RESULTS_DIR, "tx2gene.csv"))
head(tx2gene)

if (TESTING) {
  print("***Running test with Aip02 only***")
  samples <- testing_samples
} else {
  samples <- read.csv(file.path(AIPTASIA_DIR, "Samples.csv"), header=TRUE)
}
head(samples)

# Load Annotation files 
koFile <- "/work/courses/BINF6309/data_BINF6309/Module4/Annotation/ko"
pathFile <- "/work/courses/BINF6309/data_BINF6309/Module4/Annotation/path.txt"

# Load annotation file to tables 
path <- read.table(pathFile, sep="\t", header=FALSE)
ko <- read.table(koFile, sep="\t", header=FALSE)

# Set column names for each
colnames(path) <- c("ko", "pathway")
colnames(ko) <- c("pathway", "description")

# pathway identifier use only path:ko, knowledge from BINF6308 module10
path<-path %>% filter(grepl('path:ko', pathway))

# Check each tables
kable(head(path))
kable(head(ko))

files <- file.path(RESULTS_DIR, "quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = samples, 
                                design = ~ Menthol + Vibrio)

dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
  if(result != 'Intercept'){
    res <- results(dds, alpha=.05, name=result)
    dfRes <- as.data.frame(res)
    dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
    dfRes$Factor <- result
    dfRes$ko<-rownames(dfRes)
    dfAll <- rbind(dfAll, dfRes)
  }
}
# I move rbind out from the loop to avoid ko duplication.
# ex:K23034 becomes K23034 and K230341 
# dfAll <- rbind(dfAll, dfRes)
# change rowname from ko to indexes 
# dfAll <- tibble::rownames_to_column(dfAll, "ko")

# Filter the data adjusted p-value less than .05
padj05<-filter(dfAll,padj<0.05)

# Merge annotations by left join
mergeKoPath<-merge(padj05,path,by="ko",all.x = TRUE)
mergeAll<-subset(merge(mergeKoPath,ko,by="pathway",all.x = TRUE),select=c("ko","pathway","description","log2FoldChange","padj","Factor"))
# Drop nan value rows
deAnnotated<-na.omit(mergeAll)
deAnnotated
# Write table to csv file
write.csv(deAnnotated, file=file.path(RESULTS_DIR, "deAnnotated.csv"))
# end of de.R script

