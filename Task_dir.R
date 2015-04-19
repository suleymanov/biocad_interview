library(igraph)
source("graph_funcs.R")

# count <- 1500
count <- 1500
n <- 20
set.seed(2015)
gfs <- gen_seq_graphs(count, n, "dir", 0.1)

print("Run standard procedure...")
print_isomorphs(gfs)
print("Done standard procedure!")

print("Run algorithm with kernel...")
for (i in 1:(length(gfs) - 1)) {
    k_ii <- kernel1(gfs[[i]], gfs[[i]], 2, "dir")
    for (j in (i + 1):length(gfs)) {
        if (any((sort(gfs[[i]][[5]]) == sort(gfs[[j]][[5]])) == FALSE)) { next }
        k_jj <- kernel1(gfs[[j]], gfs[[j]], 2, "dir")
        if (k_ii == k_jj)
            print(c(i, j))
    }
}

gfs1 <- gen_seq_graphs(count - 50, n, "dir", 0.1)
gfs2 <- sample(gfs1, 50)
gfs <- c(gfs1, gfs2)

print("Run standard procedure...")
print_isomorphs(gfs)
print("Done standard procedure!")

print("Run algorithm with kernel...")
for (i in 1:(length(gfs) - 1)) {
    k_ii <- kernel1(gfs[[i]], gfs[[i]], 2, "dir")
    for (j in (i + 1):length(gfs)) {
        if (any((sort(gfs[[i]][[5]]) == sort(gfs[[j]][[5]])) == FALSE)) { next }
        k_jj <- kernel1(gfs[[j]], gfs[[j]], 2, "dir")
        if (k_ii == k_jj)
            print(c(i, j))
    }
}