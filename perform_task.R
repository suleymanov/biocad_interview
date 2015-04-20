library(igraph)
source("graph_funcs.R")

perform_task <- function(gfs, mode = c("dir", "undir"), f_names) {
    print("Running standard procedure...")
    print_non_isomorphs(gfs, f_names[[1]])
    print("Done with standard procedure!")
    print("")
    print("Running algorithm with kernel...")
    
    ind <- 1
    v <- list()
    for (i in 1:(length(gfs) - 1)) {
        ne_i <- ecount(gfs[[i]][[1]])
        nv_i <- vcount(gfs[[i]][[1]])
        k_ii <- kernel1(gfs[[i]], gfs[[i]], 3, mode)
        for (j in (i + 1):length(gfs)) {
            ne_j <- ecount(gfs[[j]][[1]])
            nv_j <- vcount(gfs[[j]][[1]])
            if (ne_i != ne_j) { v <- append(v, ind, i, j); ind <- ind + 1; next }
            if (nv_i != nv_j) { v <- append(v, ind, i, j); ind <- ind + 1; next }
            if (!isTRUE(all.equal(sort(gfs[[i]][[5]]), sort(gfs[[j]][[5]])))) {
                v <- append(v, ind, i, j)
                ind <- ind + 1
                next
            }
            k_jj <- kernel1(gfs[[j]], gfs[[j]], 3, mode)
            if (!isTRUE(all.equal(k_ii, k_jj))) {
                v <- append(v, ind, i, j)
                ind <- ind + 1
            }
        }
    }
    v <- as.character(v)
    write(v, file = f_names[[2]], append = F, sep = ",", ncolumns = 10)
    print("Done!")
}

append <- function(v, ind, i, j) {
    v[[ind]] <- paste("(", i, " ", j, ")", sep = "")
    return (v)
}