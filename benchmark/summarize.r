#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
infile <- args[1]

x <- read.csv(infile, header=T, sep="\t")
xs <- aggregate(x[,-1], list(x[,1]), mean)
names(xs)[1] <- "method"

write.table(xs, file=paste("mean_", infile, sep=""), append=F, quote=F, sep="\t", row.names=F, col.names=T)
