library(igraph)
# library(gputools)

# generates undirected graph with 'n' nodes
# given fraction 'prob' of all possible connections
gen_graph_undir <- function(n, prob) {
    entries <- sample(c(0, 1), sum(1:(n - 1)), replace = T, prob = c(1 - prob, prob))
    adj <- matrix(rep(0, n * n), nrow = n, ncol = n)
    ind <- 1
    for (i in 1:(n - 1)) {
        beg <- ind
        end <- ind + (n - i) - 1
        row <- entries[beg:end]
        for (j in 1:length(row)) {
            if (row[j] == 1) {
                adj[i, j + i] <- 1
                adj[j + i, i] <- 1
            }
        }
        ind <- end + 1
    }
    grf <- graph.adjacency(adj, "undirected")
    p <- matrix(rep(1 / n, n), nrow = n)
    q <- p
    P <- adj
    rs <- colSums(P)
    inds <- rs != 0
    P[, inds] <- t(P[inds, ] / rs[inds])
    fr <- 1 - sum(q[!inds])
    q[inds] <- fr * rs[inds] / sum(rs)
    res <- list(grf, adj, P, p, q)
}

# generates directed graph with 'n' nodes
# given probability 'prob' of node of having edge with any other arbitrary node
gen_graph_dir <- function(n, prob) {
    entries <- sample(c(0, 1), n * (n - 1), replace = T, prob = c(1 - prob, prob))
    adj <- matrix(rep(0, n * n), nrow = n, byrow = T)
    ind <- 1
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                adj[i, j] <- entries[ind]
                ind <- ind + 1
            }
        }
    }
    grf <- graph.adjacency(adj, "directed")
    p <- matrix(rep(1 / n, n), nrow = n)
    q <- p
    P <- adj
    rs <- rowSums(P); cs <- colSums(P)
    inds_r <- rs != 0; inds_c <- cs != 0
    P[, inds_r] <- t(P[inds_r, ] / rs[inds_r])
    # fr <- 1 - sum(q[!inds_c])
    # q[inds_c] <- fr * cs[inds_c] / sum(cs)
    fr <- 1 - sum(q[!inds_r])
    q[inds_r] <- fr * cs[inds_r] / sum(cs)
    res <- list(grf, adj, P, p, q)
}

# direct product of two graphs adjacency matrices
direct_product <- function(g1, g2, mode = c("dir", "undir")) {
    m <- kronecker(g1[[2]], g2[[2]])
    if (mode == "dir")
        g <- graph.adjacency(m, "directed")
    else
        g <- graph.adjacency(m, "undirected")
    res <- list(g, m)
}

kernel1 <- function(g1, g2, k, mode = c("dir", "undir")) {
    g <- direct_product(g1, g2, mode)
    px <- kronecker(g1[[4]], g2[[4]])
    qx <- kronecker(g1[[5]], g2[[5]])
    s <- t(qx) %*% px
    val <- diag(rep(1, dim(g[[2]])[1]))
    for (i in 1:k) {
        val <- val %*% g[[2]]
        # val <- gpuMatMult(val, g[[2]])
        s <- s + t(qx) %*% val %*% px
    }
    return (s)
}

print_k_kernels <- function(g1, g2, k, mode = c("dir", "undir")) {
    for (i in 1:k) {
        k11 <- kernel1(g1, g1, i, mode)
        k12 <- kernel1(g1, g2, i, mode)
        k22 <- kernel1(g2, g2, i, mode)
        print(c(k11, k12, k22))
    }
}

# generates sequence of 'count' graphs with 'n' nodes
gen_seq_graphs <- function(count, n, mode = c("dir", "undir"), prob) {
    grfs <- list()
    if (mode == "dir")
        for (i in 1:count)
            grfs[[i]] <- gen_graph_dir(n, prob)
    else
        for (i in 1:count)
            grfs[[i]] <- gen_graph_undir(n, prob)
    grfs
}

# prints pairs of isomorphic graphs using standard methods
print_isomorphs <- function(grfs) {
    for (i in 1:(length(grfs) - 1)) 
        for (j in (i + 1):length(grfs))
            if (graph.isomorphic(grfs[[i]][[1]], grfs[[j]][[1]]))
                print(c(i, j))
}