# Random orthogonal matrix

Generate a random orthogonal basis matrix of dimension \\((nrow x
ncol)\\) using the method in Heiberger (1978).

## Usage

``` r
randomOrthogonalMatrix(nrow, ncol, n = nrow, d = ncol, seed = NULL)
```

## Arguments

- nrow:

  the number of rows of the resulting orthogonal matrix.

- ncol:

  the number of columns of the resulting orthogonal matrix.

- n:

  deprecated. See `nrow` above.

- d:

  deprecated. See `ncol` above.

- seed:

  an optional integer argument to use in
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for
  reproducibility. By default the current seed will be used.
  Reproducibility can also be achieved by calling
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before calling this
  function.

## Details

The use of arguments `n` and `d` is deprecated and they will be removed
in the future.

## Value

An orthogonal matrix of dimension \\(nrow x ncol\\) such that each
column is orthogonal to the other and has unit lenght. Because of the
latter, it is also called orthonormal.

## See also

[`coordProj`](https://mclust-org.github.io/mclust/reference/coordProj.md)

## References

Heiberger R. (1978) Generation of random orthogonal matrices. *Journal
of the Royal Statistical Society. Series C (Applied Statistics)*, 27(2),
199-206.

## Examples

``` r
B <- randomOrthogonalMatrix(10,3)
zapsmall(crossprod(B))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
```
