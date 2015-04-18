# generates undirected graph with 'n' nodes
# given fraction 'prob' of all possible connections
gen_graph_undir <- function(n, prob) {
    graph_sz <- sum(1:n)
    entries <- sample(c(0, 1), sum(1:(n - 1)), replace = T, prob = c(1 - prob, prob))
    m <- matrix(rep(0, n * n), nrow = n, ncol = n)
    ind <- 1
    for (i in 1:(n - 1)) {
        beg <- ind
        end <- ind + (n - i) - 1
        row <- entries[beg:end]
        for (j in 1:length(row)) {
            if (row[j] == 1) {
                m[i, j + i] <- 1
                m[j + i, i] <- 1
            }
        }
        ind <- end + 1
    }
    g <- graph.adjacency(m, "undirected")
    res <- list(m, g)
}

num_e <- function(g, mode = c("dir", "undir")) {
    if (mode == "dir") return (sum(g[[1]][[1]]))
    else return (sum(g[[1]]) / 2.0)
}

# generates directed graph with 'n' nodes
# given probability 'prob' of node of having edge with any other arbitrary node
gen_graph_dir <- function(n, prob) {
    entries <- sample(c(0, 1), n * (n - 1), replace = T, prob = c(1 - prob, prob))
    m <- matrix(rep(0, n * n), nrow = n, byrow = T)
    ind <- 1
    for (i in 1:n) {
        for (j in 1:n) {
            if (i != j) {
                m[i, j] <- entries[ind]
                ind <- ind + 1
            }
        }
    }
    g <- graph.adjacency(m, "directed")
    res <- list(m, g)
}

# direct product of two graphs adjacency matrices
direct_product <- function(g1, g2, mode = c("dir", "undir")) {
    m <- kronecker(g1[[1]], g2[[1]])
    if (mode == "dir")
        g <- graph.adjacency(m, "directed")
    else
        g <- graph.adjacency(m, "undirected")
    res <- list(m, g)
}

kernel1 <- function(Wx, k) {
    s <- t(qx) %*% px
    val <- diag(rep(1, dim(Wx)[1]))
    for (i in 1:k) {
        val <- val %*% Wx
        s <- s + t(qx) %*% val %*% px
    }
    s
}

print_k_kernels <- function(g1, g2, k, mode = c("dir", "undir")) {
    g11 <- direct_product(g1, g1, mode)
    g12 <- direct_product(g1, g2, mode)
    g22 <- direct_product(g2, g2, mode)
    for (i in 1:k) { print(c(kernel1(g11[[1]], i), kernel1(g12[[1]], i), kernel1(g22[[1]], i))) }
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
    for (i in 1:(length(grfs) - 1)) {
        for (j in (i + 1):length(grfs)) {
            if (graph.isomorphic(grfs[[i]][[2]], grfs[[j]][[2]])) {
                print(c(i, j))
            }
        }
    }
}