# Adjusted Rand Index

Computes the adjusted Rand index comparing two classifications.

## Usage

``` r
adjustedRandIndex(x, y)
```

## Arguments

- x:

  A numeric or character vector of class labels.

- y:

  A numeric or character vector of class labels. The length of `y`
  should be the same as that of `x`.

## Value

The adjusted Rand index comparing the two partitions (a scalar). This
index has zero expected value in the case of random partition, and it is
bounded above by 1 in the case of perfect agreement between two
partitions.

## References

L. Hubert and P. Arabie (1985) Comparing Partitions, *Journal of the
Classification*, 2, pp. 193-218.

## See also

[`classError`](https://mclust-org.github.io/mclust/reference/classError.md),
[`mapClass`](https://mclust-org.github.io/mclust/reference/mapClass.md),
[`table`](https://rdrr.io/r/base/table.html)

## Examples

``` r
a <- rep(1:3, 3)
a
#> [1] 1 2 3 1 2 3 1 2 3
b <- rep(c("A", "B", "C"), 3)
b
#> [1] "A" "B" "C" "A" "B" "C" "A" "B" "C"
adjustedRandIndex(a, b)
#> [1] 1

a <- sample(1:3, 9, replace = TRUE)
a
#> [1] 3 2 3 3 2 3 1 2 2
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
#> [1] "A" "A" "B" "A" "B" "C" "C" "B" "B"
adjustedRandIndex(a, b)
#> [1] 0.08695652

a <- rep(1:3, 4)
a
#>  [1] 1 2 3 1 2 3 1 2 3 1 2 3
b <- rep(c("A", "B", "C", "D"), 3)
b
#>  [1] "A" "B" "C" "D" "A" "B" "C" "D" "A" "B" "C" "D"
adjustedRandIndex(a, b)
#> [1] -0.2790698

irisHCvvv <- hc(modelName = "VVV", data = iris[,-5])
cl3 <- hclass(irisHCvvv, 3)
adjustedRandIndex(cl3,iris[,5])
#> [1] 0.7591987

irisBIC <- mclustBIC(iris[,-5])
adjustedRandIndex(summary(irisBIC,iris[,-5])$classification,iris[,5])
#> [1] 0.5681159
adjustedRandIndex(summary(irisBIC,iris[,-5],G=3)$classification,iris[,5])
#> [1] 0.9038742
```
