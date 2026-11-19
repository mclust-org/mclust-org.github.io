# Convert mixture component covariances to decomposition form.

Converts a set of covariance matrices from representation as a 3-D array
to a parameterization by eigenvalue decomposition.

## Usage

``` r
sigma2decomp(sigma, G = NULL, tol = sqrt(.Machine$double.eps), ...)
```

## Arguments

- sigma:

  Either a 3-D array whose \[,,k\]th component is the covariance matrix
  for the kth component in an MVN mixture model, or a single covariance
  matrix in the case that all components have the same covariance.

- G:

  The number of components in the mixture. When `sigma` is a 3-D array,
  the number of components can be inferred from its dimensions.

- tol:

  Tolerance for determining whether or not the covariances have equal
  volume, shape, and or orientation. The default is the square root of
  the relative machine precision, `sqrt(.Machine$double.eps)`, which is
  about `1.e-8`.

- ...:

  Catches unused arguments from an indirect or list call via `do.call`.

## Value

The covariance matrices for the mixture components in decomposition
form, including the following components:

- modelName:

  A character string indicating the infered model. The help file for
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)
  describes the available models.

- d:

  The dimension of the data.

- G:

  The number of components in the mixture model.

- scale:

  Either a *G*-vector giving the scale of the covariance (the *d*th root
  of its determinant) for each component in the mixture model, or a
  single numeric value if the scale is the same for each component.

- shape:

  Either a *G* by *d* matrix in which the *k*th column is the shape of
  the covariance matrix (normalized to have determinant 1) for the *k*th
  component, or a *d*-vector giving a common shape for all components.

- orientation:

  Either a *d* by *d* by *G* array whose `[,,k]`th entry is the
  orthonomal matrix whose columns are the eigenvectors of the covariance
  matrix of the *k*th component, or a *d* by *d* orthonormal matrix if
  the mixture components have a common orientation. The `orientation`
  component of `decomp` can be omitted in spherical and diagonal models,
  for which the principal components are parallel to the coordinate axes
  so that the orientation matrix is the identity.

## See also

[`decomp2sigma`](https://mclust-org.github.io/mclust/reference/decomp2sigma.md)

## Examples

``` r
meEst <- meEEE(iris[,-5], unmap(iris[,5])) 
names(meEst$parameters$variance)
#> [1] "modelName" "d"         "G"         "sigma"     "Sigma"     "cholSigma"
meEst$parameters$variance$Sigma
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.26387916  0.08987702   0.16956251  0.03932391
#> Sepal.Width    0.08987702  0.11197340   0.05118186  0.03002534
#> Petal.Length   0.16956251  0.05118186   0.18637706  0.04196406
#> Petal.Width    0.03932391  0.03002534   0.04196406  0.03974960

sigma2decomp(meEst$parameters$variance$Sigma, G = length(unique(iris[,5])))
#> $sigma
#> , , 1
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.26387916 0.08987702 0.16956251 0.03932391
#> [2,] 0.08987702 0.11197340 0.05118186 0.03002534
#> [3,] 0.16956251 0.05118186 0.18637706 0.04196406
#> [4,] 0.03932391 0.03002534 0.04196406 0.03974960
#> 
#> , , 2
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.26387916 0.08987702 0.16956251 0.03932391
#> [2,] 0.08987702 0.11197340 0.05118186 0.03002534
#> [3,] 0.16956251 0.05118186 0.18637706 0.04196406
#> [4,] 0.03932391 0.03002534 0.04196406 0.03974960
#> 
#> , , 3
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.26387916 0.08987702 0.16956251 0.03932391
#> [2,] 0.08987702 0.11197340 0.05118186 0.03002534
#> [3,] 0.16956251 0.05118186 0.18637706 0.04196406
#> [4,] 0.03932391 0.03002534 0.04196406 0.03974960
#> 
#> 
#> $d
#> [1] 4
#> 
#> $modelName
#> [1] "EEE"
#> 
#> $G
#> [1] 3
#> 
#> $scale
#> [1] 0.08114033
#> 
#> $shape
#> [1] 5.4441760 1.0696919 0.6344826 0.2706384
#> 
#> $orientation
#>            [,1]        [,2]       [,3]       [,4]
#> [1,] -0.7401980  0.01338854  0.6352049  0.2200963
#> [2,] -0.3056652 -0.87119601 -0.2316293 -0.3064868
#> [3,] -0.5783262  0.47573614 -0.5625433 -0.3503698
#> [4,] -0.1556100 -0.12047151 -0.4758232  0.8572423
#> 
```
