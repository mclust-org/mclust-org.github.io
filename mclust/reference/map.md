# Classification given Probabilities

Converts a matrix in which each row sums to 1 to an integer vector
specifying for each row the column index of the maximum.

## Usage

``` r
map(z, warn = mclust.options("warn"), ...)
```

## Arguments

- z:

  A matrix (for example a matrix of conditional probabilities in which
  each row sums to 1 as produced by the E-step of the EM algorithm).

- warn:

  A logical variable indicating whether or not a warning should be
  issued when there are some columns of `z` for which no row attains a
  maximum.

- ... :

  Provided to allow lists with elements other than the arguments can be
  passed in indirect or list calls with `do.call`.

## Value

A integer vector with one entry for each row of z, in which the *i*-th
value is the column index at which the *i*-th row of `z` attains a
maximum.

## See also

[`unmap`](https://mclust-org.github.io/mclust/reference/unmap.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`em`](https://mclust-org.github.io/mclust/reference/em.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md).

## Examples

``` r
emEst <- me(modelName = "VVV", data = iris[,-5], z = unmap(iris[,5]))

map(emEst$z)
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 3 2 3 2
#>  [75] 2 2 2 3 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3
#> [112] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
#> [149] 3 3
```
