# Softmax function

Efficient implementation (via Fortran) of the softmax (aka multinomial
logistic) function converting a set of numerical values to probabilities
summing to 1.

## Usage

``` r
softmax(x, v = NULL)
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
x_k\]\\) (with \\(i=1,\dots,n\\)), the softmax function calculates \$\$
\text{softmax}(x\_{\[i\]})\_j = \dfrac{\exp(x_j + v_j)}{\sum\_{l=1}^k
\exp(x_l + v_l)} \qquad \text{for } j = 1,\dots,k \$\$

## Value

Returns a matrix of the same dimension as `x` with values in the range
\\((0,1)\\) that sum to 1 along the rows.

## Author

Luca Scrucca

## See also

[`logsumexp`](https://mclust-org.github.io/mclust/reference/logsumexp.md)

## References

Blanchard P., Higham D. J., Higham N. J. (2021). Accurately computing
the log-sum-exp and softmax functions. *IMA Journal of Numerical
Analysis*, 41/4:2311–2330.
[doi:10.1093/imanum/draa038](https://doi.org/10.1093/imanum/draa038)

## Examples

``` r
x = matrix(rnorm(15), 5, 3)
v = log(c(0.5, 0.3, 0.2))
(z = softmax(x, v))
#>           [,1]      [,2]       [,3]
#> [1,] 0.2554222 0.7150265 0.02955130
#> [2,] 0.1034917 0.8674056 0.02910264
#> [3,] 0.3018665 0.3382748 0.35985864
#> [4,] 0.7205435 0.2329123 0.04654415
#> [5,] 0.6109547 0.2555151 0.13353026
rowSums(z)
#> [1] 1 1 1 1 1
```
