# Weighted means, covariance and scattering matrices conditioning on a weighted matrix

Compute efficiently (via Fortran code) the means, covariance and
scattering matrices conditioning on a weighted or indicator matrix

## Usage

``` r
covw(X, Z, normalize = TRUE)
```

## Arguments

- X:

  A \\((n x p)\\) data matrix, with \\(n\\) observations on \\(p\\)
  variables.

- Z:

  A \\((n x G)\\) matrix of weights, with \\(G\\) number of groups.

- normalize:

  A logical indicating if rows of `Z` should be normalized to sum to
  one.

## Value

A list with the following components:

- mean:

  A \\((p x G)\\) matrix of weighted means.

- S:

  A \\((p x p x G)\\) array of weighted covariance matrices.

- W:

  A \\((p x p x G)\\) array of weighted scattering matrices.

## Author

M. Fop and L. Scrucca

## Examples

``` r
# Z as an indicator matrix
X <- iris[,1:4]
Z <- unmap(iris$Species)
str(covw(X, Z))
#> List of 3
#>  $ mean: num [1:4, 1:3] 5.006 3.428 1.462 0.246 5.936 ...
#>  $ S   : num [1:4, 1:4, 1:3] 0.1218 0.0972 0.016 0.0101 0.0972 ...
#>  $ W   : num [1:4, 1:4, 1:3] 6.088 4.862 0.801 0.506 4.862 ...
# Z as a matrix of weights
mod <- Mclust(X, G = 3, modelNames = "VVV")
str(covw(X, mod$z))
#> List of 3
#>  $ mean: num [1:4, 1:3] 5.006 3.428 1.462 0.246 5.915 ...
#>  $ S   : num [1:4, 1:4, 1:3] 0.1218 0.0972 0.016 0.0101 0.0972 ...
#>  $ W   : num [1:4, 1:4, 1:3] 6.088 4.862 0.801 0.506 4.862 ...
```
