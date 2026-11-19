# Indicator Variables given Classification

Converts a classification into a matrix of indicator variables.

## Usage

``` r
unmap(classification, groups=NULL, noise=NULL, ...)
```

## Arguments

- classification:

  A numeric or character vector. Typically the distinct entries of this
  vector would represent a classification of observations in a data set.

- groups:

  A numeric or character vector indicating the groups from which
  `classification` is drawn. If not supplied, the default is to assumed
  to be the unique entries of classification.

- noise:

  A single numeric or character value used to indicate the value of
  `groups` corresponding to noise.

- ... :

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

An *n* by *m* matrix of *(0,1)* indicator variables, where *n* is the
length of `classification` and *m* is the number of unique values or
symbols in `classification`. Columns are labeled by the unique values in
`classification`, and the `[i,j]`th entry is *1* if `classification[i]`
is the *j*th unique value or symbol in sorted order `classification`. If
a `noise` value of symbol is designated, the corresponding indicator
variables are relocated to the last column of the matrix.

## See also

[`map`](https://mclust-org.github.io/mclust/reference/map.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md)

## Examples

``` r
z <- unmap(iris[,5])
z[1:5, ]
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    1    0    0
#> [3,]    1    0    0
#> [4,]    1    0    0
#> [5,]    1    0    0
  
emEst <- me(modelName = "VVV", data = iris[,-5], z = z)
emEst$z[1:5,]
#>      [,1]         [,2]         [,3]
#> [1,]    1 1.340380e-44 1.861339e-34
#> [2,]    1 2.201405e-31 6.676298e-28
#> [3,]    1 1.896748e-36 1.102178e-29
#> [4,]    1 3.488647e-32 6.409600e-26
#> [5,]    1 4.393475e-47 7.745885e-35
  
map(emEst$z)
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 2 3 2 3 2
#>  [75] 2 2 2 3 2 2 2 2 2 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3
#> [112] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
#> [149] 3 3
```
