# Log sum of exponentials

Efficient implementation (via Fortran) of the log-sum-exp function.

## Usage

``` r
logsumexp(x, v = NULL)
```

## Arguments

- x:

  a matrix of dimension \\(n \times k\\) of numerical values. If a
  vector is provided, it is converted to a single-row matrix.

- v:

  an optional vector of length \\(k\\) of numerical values to be added
  to each row of `x` matrix. If not provided, a vector of zeros is used.

## Details

Given the matrix `x`, for each row \\(x\_{\[i\]} = \[x_1, \dots,
x_k\]\\) (with \\(i=1,\dots,n\\)), the log-sum-exp (LSE) function
calculates \$\$ \text{LSE}(x\_{\[i\]}) = \log \sum\_{j=1}^k \exp(x_j +
v_j) = m + \log \sum\_{j=1}^k \exp(x_j + v_j - m) \$\$ where \\(m =
\max(x_1+v_1, \dots, x_k+v_k)\\).

## Value

Returns a vector of values of length equal to the number of rows of `x`.

## Author

Luca Scrucca

## See also

[`softmax`](https://mclust-org.github.io/mclust/reference/softmax.md)

## References

Blanchard P., Higham D. J., Higham N. J. (2021). Accurately computing
the log-sum-exp and softmax functions. *IMA Journal of Numerical
Analysis*, 41/4:2311–2330.
[doi:10.1093/imanum/draa038](https://doi.org/10.1093/imanum/draa038)

## Examples

``` r
x = matrix(rnorm(15), 5, 3)
v = log(c(0.5, 0.3, 0.2))
logsumexp(x, v)
#> [1] -0.6711947 -0.6393459  0.5306414 -0.1252858  0.8242612
```
