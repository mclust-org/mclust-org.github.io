# Classifies Data According to Unique Observations

Gives a one-to-one mapping from unique observations to rows of a data
matrix.

## Usage

``` r
partuniq(x)
```

## Arguments

- x:

  Matrix of observations.

## Value

A vector of length `nrow(x)` with integer entries. An observation `k` is
assigned an integer `i` whenever observation `i` is the first row of `x`
that is identical to observation `k` (note that `i <= k`).

## See also

[`partconv`](https://mclust-org.github.io/mclust/reference/partconv.md)

## Examples

``` r
set.seed(0)

mat <- data.frame(lets = sample(LETTERS[1:2],9,TRUE), nums = sample(1:2,9,TRUE))
mat
#>   lets nums
#> 1    B    2
#> 2    A    2
#> 3    B    1
#> 4    A    1
#> 5    A    1
#> 6    B    1
#> 7    A    1
#> 8    A    2
#> 9    A    2

ans <- partuniq(mat)
ans
#> [1] 1 2 3 4 4 3 4 2 2

partconv(ans,consec=TRUE)
#> [1] 1 2 3 4 4 3 4 2 2
```
