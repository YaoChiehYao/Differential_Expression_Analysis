##  BINF 6309 Assignment 4 Differential expression analysis
Author: Yao Chieh Yao


## Description
In assignment 4, this experiment applies Aiptasia Miseq data for 
differential expression analysis in two main steps. First, build
indexes and align the reads of multiple FASTQ files by Salmon 
to estimate the relative abundance within the referenced Aiptasia 
genome. Second, apply tximport library in R script to import the 
Salmon abundance estimates and DESeq2 package to perform statistical 
tests to identify deferentially expressed genes. 


## Getting Started
* Hi, this is the documentation for assignment four of the bio-computational
  method course, BINF6309.
* Please download the file in the GitHub link below and follow the
  executing section's guidelines in the Rmarkdown file. The working environment 
  recommends the command line prompt and Rstudio. 
 
  Here is the link to the file: 
```
https://github.com/NU-Bioinformatics/module04-YaoChiehYao.git
```


## Method and Result
Please read my Rmarkdown "methodsResults.Rmd" in the GitHub link above 
chunk by chuck to understand the RNA-seq analysis pipeline and the 
comparative methods used in this assignment.   


## References
Love, Michael I., Wolfgang Huber, and Simon Anders. 2014. “Moderated Estimation of 
Fold Change and Dispersion for RNA-Seq Data with DESeq2.” *Genome Biol* 15 (12): 550–50.

Patro, Rob, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and Carl Kingsford. 
2017. “Salmon Provides Fast and Bias-Aware Quantification of Transcript Expression.” 
*Nat Methods* 14 (4): 417–19.

Soneson, Charlotte, Michael I. Love, and Mark D. Robinson. 2016. “Differential Analyses 
for RNA-Seq: Transcript-Level Estimates Improve Gene-Level Inferences.” *F1000Res* 4 
(February): 1521–1.
