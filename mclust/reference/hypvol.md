# Approximate Hypervolume for Multivariate Data

Computes a simple approximation to the hypervolume of a multivariate
data set.

## Usage

``` r
hypvol(data, reciprocal=FALSE)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- reciprocal:

  A logical variable indicating whether or not the reciprocal
  hypervolume is desired rather than the hypervolume itself. The default
  is to return the hypervolume.

## Value

Returns the minimum of the hypervolume computed from simple variable
bounds and that computed from variable bounds of the principal component
scores. Used for the default hypervolume parameter for the noise
component when observations are designated as noise in `Mclust` and
`mclustBIC`.

## References

A. Dasgupta and A. E. Raftery (1998). Detecting features in spatial
point processes with clutter via model-based clustering. *Journal of the
American Statistical Association 93:294-302*.

C. Fraley and A.E. Raftery (1998). *Computer Journal 41:578-588*.

C. Fraley and A. E. Raftery (2002). Model-based clustering, discriminant
analysis, and density estimation. *Journal of the American Statistical
Association 97:611-631*.

## See also

[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)

## Examples

``` r
hypvol(iris[,-5])
#> [1] 27.2412
```
