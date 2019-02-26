#!/usr/bin/env Rscript
#Author: Pulkit Anupam Srivastava
#Version: 1.0
#Last Modified: 26 Feb, 2019
#Description: The script generate feature-sets for given protein sequence(s).
#Input: Path to fasta file. Name of fasta file.
#Output: Generates files containing descriptor value for all sequences in fasta file.
library("protr")
library("foreign")
#
args = commandArgs(TRUE)
path2dir <- args[1]
FastaFile <- args[2]
#
seq_class1 = readFASTA(FastaFile)
seq_class1 = seq_class1[(sapply(seq_class1, protcheck))]
#List of feature-sets
func2perform=c('extractMoreauBroto','extractTC','extractDC',
               'extractMoran','extractGeary','extractCTriad')
#Iteration over list of feature-sets
for (f in func2perform)
{
  class1 = t(sapply(seq_class1, f))
  #Create a folder with feature-set name
  dir.create(file.path(path2dir, f), showWarnings = FALSE)
  #
  outputfile_csv=file.path(path2dir,f,paste(f,".csv",sep=""))
  write.csv(class1,outputfile_csv,sep=",",row.names=TRUE)
  #Writes csv file having descriptor values for each protein sequence
  data_class1 = read.csv(outputfile_csv,header=TRUE)
  colnames(data_class1)[1]<-"ID"
  data_class1$label <- "?"
  mydata=rbind(data_class1)
  write.csv(mydata,outputfile_csv,sep=",",row.names=TRUE)
}
