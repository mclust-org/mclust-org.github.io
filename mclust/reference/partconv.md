# Numeric Encoding of a Partitioning

Converts a vector interpreted as a classification or partitioning into a
numeric vector.

## Usage

``` r
partconv(x, consec=TRUE)
```

## Arguments

- x:

  A vector interpreted as a classification or partitioning.

- consec:

  Logical value indicating whether or not to consecutive class numbers
  should be used .

## Value

Numeric encoding of `x`. When `consec = TRUE`, the distinct values in
`x` are numbered by the order in which they appear. When
`consec = FALSE`, each distinct value in `x` is numbered by the index
corresponding to its first appearance in `x`.

## See also

[`partuniq`](https://mclust-org.github.io/mclust/reference/partuniq.md)

## Examples

``` r
partconv(iris[,5])
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [75] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3
#> [112] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
#> [149] 3 3

set.seed(0)
cl <- sample(LETTERS[1:9], 25, replace=TRUE)
partconv(cl, consec=FALSE)
#>  [1]  1  2  3  4  5  3  5  8  4 10 10 12  3  1 10 10  1  1 10 10  5  1  4  2  8
partconv(cl, consec=TRUE)
#>  [1] 1 2 3 4 5 3 5 6 4 7 7 8 3 1 7 7 1 1 7 7 5 1 4 2 6
```
