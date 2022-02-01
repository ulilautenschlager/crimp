#!/usr/bin/env Rscript

library(pophelper)
x <- readQ("input_pophelper.txt")

qmat <- lapply(x, as.matrix)
p <- aperm(simplify2array(qmat), c(3, 1, 2))
perm <- label.switching::stephens(p)

write.table(perm$permutations, file="permutations.csv", quote=F, row.names=F, col.names=F)
