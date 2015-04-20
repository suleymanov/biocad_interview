library(igraph)
source("graph_funcs.R")

perform_task <- function(gfs, mode = c("dir", "undir"), f_names) {
    print("Running standard procedure...")
    v_ref <- print_non_isomorphs(gfs, f_names[[1]])
    print("Done with standard procedure!")
    print("")
    print("Running algorithm with kernel...")
    
    ind <- 1
    v <- list()
    v_ch <- list()
    for (i in 1:(length(gfs) - 1)) {
        ne_i <- ecount(gfs[[i]][[1]])
        nv_i <- vcount(gfs[[i]][[1]])
        k_ii <- kernel2(gfs[[i]], 3, mode)
        for (j in (i + 1):length(gfs)) {
            ne_j <- ecount(gfs[[j]][[1]])
            nv_j <- vcount(gfs[[j]][[1]])
            if (ne_i != ne_j || nv_i != nv_j) {
                res <- append(ind, i, j)
                v[[ind]] <- res[[1]]
                v_ch[[ind]] <- res[[2]]
                ind <- ind + 1
                next
            }
            if (!isTRUE(all.equal(sort(gfs[[i]][[5]]), sort(gfs[[j]][[5]])))) {
                res <- append(ind, i, j)
                v[[ind]] <- res[[1]]
                v_ch[[ind]] <- res[[2]]
                ind <- ind + 1
                next
            }
            k_jj <- kernel2(gfs[[j]], 3, mode)
            if (!isTRUE(all.equal(k_ii, k_jj))) {
                res <- append(ind, i, j)
                v[[ind]] <- res[[1]]
                v_ch[[ind]] <- res[[2]]
                ind <- ind + 1
            }
        }
    }
    v_ch <- as.character(v_ch)
    write(v_ch, file = f_names[[2]], append = F, sep = ",", ncolumns = 10)
    print("Done!")
    
    print("Diagnosis: ")
    print(length(v) == length(v_ref))
    print(all.equal(v, v_ref))
}

append <- function(ind, i, j) {
    v <- c(i, j)
    v_ch <- paste("(", i, " ", j, ")", sep = "")
    res <- list(v, v_ch)
}

