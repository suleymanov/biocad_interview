# Some playing with graphs #2
Rail Suleymanov  
Friday, April 17, 2015  

This is a demo document on graphs - generations and some calculations (case of directed graphs).

We begin with loading graph library and some useful functions:


```r
library(igraph)
source("graph_funcs.R")
```

Then we generate 20 graphs, each consisting of 4 nodes, with probability of 0.5 of having and edge between any arbitrary edges. We will also define uniform initial and stopping probabilities for all graphs and calculate direct product of those probabilities:


```r
set.seed(2015)
gfs <- gen_seq_graphs(20, 4, "dir", 0.5)
p <- matrix(rep(0.25, 4), nrow = 4)
q <- p
px <- kronecker(p, p)
qx <- kronecker(q, q)
```

Now we try to find all pairs of isomorphs using standard methods:


```r
print_isomorphs(gfs)
```

```
## [1] 14 19
```

Plot the single pair of isomorphs:


```r
par(mfrow = c(1, 2))
plot(gfs[[14]][[2]])
plot(gfs[[19]][[2]])
```

![](Demo2_files/figure-html/unnamed-chunk-4-1.png) 

As in previous demo, we'll compute direct product and kernel value on graph #14 with itself and with graph #19 (its isomorph):


```r
g14_14 <- direct_product(gfs[[14]], gfs[[14]], "dir")
g14_19 <- direct_product(gfs[[14]], gfs[[19]], "dir")
kernel1(g14_14[[1]], 1)
```

```
##        [,1]
## [1,] 0.3125
```

```r
kernel1(g14_19[[1]], 1)
```

```
##        [,1]
## [1,] 0.3125
```

```r
kernel1(g14_14[[1]], 5)
```

```
##          [,1]
## [1,] 72.60156
```

```r
kernel1(g14_19[[1]], 5)
```

```
##          [,1]
## [1,] 72.60156
```

Let's try same approach as in Demo1 to find isomorphisms between graphs:


```r
for (i in 1:(length(gfs) - 1)) {
    g_ii <- direct_product(gfs[[i]], gfs[[i]], "dir")
    rel <- kernel1(g_ii[[1]], 3)
    for (j in (i + 1):length(gfs)) {
        g_ij <- direct_product(gfs[[i]], gfs[[j]], "dir")
        val <- kernel1(g_ij[[1]], 3)
        if (val == rel) {
            print(c(i, j))
        }
    }
}
```

```
## [1]  1 16
## [1] 10 20
## [1] 14 19
```
