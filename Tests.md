# Some tests
Rail Suleymanov  
Monday, April 20, 2015  

This document runs several tests for algorithm that finds indices of non-isomorphic graphs. Each test result will be written to 2 files - one for result, and one for reference solution. Two criteria are checked: lengths of sets of graphs indices and graph indices themselves. Double 'TRUE' declares successful test checking.

Source code files are loaded first and random seed is set:


```r
source("perform_task.R")
set.seed(2015)
```

1) First test: set of 20 undirected graphs, each of 4 vertices:


```r
gfs1 <- gen_seq_graphs(20, 4, "undir", 0.5)
perform_task(gfs1, "undir", c("20_4_undir (ref).csv", "20_4_undir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

2) Second test: set of 20 directed graphs, each of 4 vertices:


```r
gfs2 <- gen_seq_graphs(20, 4, "dir", 0.5)
perform_task(gfs2, "dir", c("20_4_dir (ref).csv", "20_4_dir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

3) Third test: set of 100 undirected graphs, each of 10 vertices:


```r
gfs3 <- gen_seq_graphs(100, 10, "undir", 0.5)
perform_task(gfs3, "undir", c("100_10_undir (ref).csv", "100_10_undir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

4) Third test: set of 100 directed graphs, each of 10 vertices:


```r
gfs4 <- gen_seq_graphs(100, 10, "dir", 0.5)
perform_task(gfs4, "dir", c("100_10_dir (ref).csv", "100_10_dir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

5) Test: set of 250 undirected graphs, each of 15 vertices:


```r
gfs5 <- gen_seq_graphs(250, 15, "undir", 0.5)
perform_task(gfs5, "undir", c("250_15_undir (ref).csv", "250_15_undir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

6) Test: set of 250 directed graphs, each of 15 vertices:


```r
gfs6 <- gen_seq_graphs(250, 15, "dir", 0.5)
perform_task(gfs6, "dir", c("250_15_dir (ref).csv", "250_15_dir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

7) Test: set of 300 undirected graphs, each of 15 vertices, with lower connection probability:


```r
gfs7 <- gen_seq_graphs(300, 15, "undir", 0.2)
perform_task(gfs7, "undir", c("300_15_undir (ref).csv", "300_15_undir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```

8) Test: set of 300 directed graphs, each of 15 vertices, with lower connection probability:


```r
gfs8 <- gen_seq_graphs(300, 15, "dir", 0.2)
perform_task(gfs8, "dir", c("300_15_dir (ref).csv", "300_15_dir.csv"))
```

```
## [1] "Running standard procedure..."
## [1] "Done with standard procedure!"
## [1] ""
## [1] "Running algorithm with kernel..."
## [1] "Done!"
## [1] "Diagnosis: "
## [1] TRUE
## [1] TRUE
```
