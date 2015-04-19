library(igraph)
source("graph_funcs2.R")

# count <- 1500
count <- 50
n <- 20
set.seed(2015)
gfs <- gen_seq_graphs(count, n, "undir", 0.5)

print("Run standard procedure...")
print("Done standard procedure!")
print_isomorphs(gfs)

print("Run algorithm with kernel...")
for (i in 1:(length(gfs) - 1)) {
    k_ii <- kernel1(gfs[[i]], gfs[[i]], 1, "undir")
    for (j in (i + 1):length(gfs)) {
        k_jj <- kernel1(gfs[[j]], gfs[[j]], 1, "undir")
        if (k_ii == k_jj)
            print(c(i, j))
    }
}

gfs1 <- gen_seq_graphs(count - 15, n, "undir", 0.5)
gfs2 <- sample(gfs1, 15)
gfs <- c(gfs1, gfs2)

print("Run standard procedure...")
print_isomorphs(gfs)
print("Done standard procedure!")

print("Run algorithm with kernel...")
for (i in 1:(length(gfs) - 1)) {
    k_ii <- kernel1(gfs[[i]], gfs[[i]], 1, "undir")
    for (j in (i + 1):length(gfs)) {
        k_jj <- kernel1(gfs[[j]], gfs[[j]], 1, "undir")
        if (k_ii == k_jj)
            print(c(i, j))
    }
}