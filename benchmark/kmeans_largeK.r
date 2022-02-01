#!/usr/bin/env Rscript

set.seed(12)
n_runs <- 100
n_clusters <- 100
n_items <- 1000  # 10 per cluster

x <- runif(n_clusters)
y <- runif(n_clusters)

# 20 points per cluster
x <- x + rnorm(n_items, sd=0.01)
y <- y + rnorm(n_items, sd=0.01)

for (run in 1:n_runs) {
    cl <- kmeans(cbind(x,y), centers=n_clusters)

    m <- matrix(0, nrow=n_items, ncol=n_clusters)
    for (i in 1:n_items) {
        m[i, cl$cluster[i]] <- 1
    }
    write.table(m, file="largeK.txt", append=ifelse(run==1, F, T), quote=F, sep=" ", row.names=F, col.names=F)
    cat("\n", file="largeK.txt", append=T)
}
