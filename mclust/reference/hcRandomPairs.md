# Random hierarchical structure

Create a hierarchical structure using a random hierarchical partition of
the data.

## Usage

``` r
hcRandomPairs(data, seed = NULL, ...)
```

## Arguments

- data:

  A numeric matrix or data frame of observations. If a matrix or data
  frame, rows correspond to observations and columns correspond to
  variables.

- seed:

  Optional single value, interpreted as an integer, specifying the seed
  for random partition.

- ...:

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A numeric two-column matrix in which the *i*th row gives the minimum
index for observations in each of the two clusters merged at the *i*th
stage of a random agglomerative hierarchical clustering.

## See also

[`hc`](https://mclust-org.github.io/mclust/reference/hc.md),
[`hclass`](https://mclust-org.github.io/mclust/reference/hclass.md)
[`hcVVV`](https://mclust-org.github.io/mclust/reference/hcE.md)

## Examples

``` r
data <- iris[,1:4]
randPairs <- hcRandomPairs(data)
str(randPairs)
#>  int [1:2, 1:149] 75 106 20 61 52 104 54 142 73 78 ...
#>  - attr(*, "initialPartition")= int [1:150] 1 2 3 4 5 6 7 8 9 10 ...
#>  - attr(*, "dimensions")= num [1:2] 150 2
# start model-based clustering from a random partition
mod <- Mclust(data, initialization = list(hcPairs = randPairs))
summary(mod)
#> ---------------------------------------------------- 
#> Gaussian finite mixture model fitted by EM algorithm 
#> ---------------------------------------------------- 
#> 
#> Mclust VVV (ellipsoidal, varying volume, shape, and orientation) model with 2
#> components: 
#> 
#>  log-likelihood   n df       BIC       ICL
#>       -214.3547 150 29 -574.0178 -574.0191
#> 
#> Clustering table:
#>   1   2 
#> 100  50 
```
