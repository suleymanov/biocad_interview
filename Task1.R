library(igraph)
source("graph_funcs.R")

# count <- 1500
count <- 100
n <- 20
set.seed(2015)
gfs <- gen_seq_graphs(count, n, "undir", 0.5) 
p <- matrix(rep(0.25, n), nrow = n) # vector of uniform starting probabilities
q <- p                              # vector of uniform stopping probabilities
px <- kronecker(p, p)
qx <- kronecker(q, q)

print("Run standard procedure...")
print("Done standard procedure!")
print_isomorphs(gfs)

print("Run algorithm with kernel...")
for (i in 1:(length(gfs) - 1)) {
    g_ii <- direct_product(gfs[[i]], gfs[[i]], "undir")
    rel <- kernel1(g_ii[[1]], 5)
    num_i <- num_e(gfs[[i]], "undir")
    for (j in (i + 1):length(gfs)) {
        num_j <- num_e(gfs[[j]], "undir")
        if (num_i == num_j) {
            g_ij <- direct_product(gfs[[i]], gfs[[j]], "undir")
            val <- kernel1(g_ij[[1]], 5)
            if (rel == val)
                print(c(i, j))
        }
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
    g_ii <- direct_product(gfs[[i]], gfs[[i]], "undir")
    rel <- kernel1(g_ii[[1]], 5)
    num_i <- num_e(gfs[[i]], "undir")
    for (j in (i + 1):length(gfs)) {
        num_j <- num_e(gfs[[j]], "undir")
        if (num_i == num_j) {
            g_ij <- direct_product(gfs[[i]], gfs[[j]], "undir")
            val <- kernel1(g_ij[[1]], 5)
            if (rel == val)
                print(c(i, j))
        }
    }
}