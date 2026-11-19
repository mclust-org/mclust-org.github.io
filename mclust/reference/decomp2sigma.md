# Convert mixture component covariances to matrix form

Converts covariances from a parameterization by eigenvalue decomposition
or cholesky factorization to representation as a 3-D array.

## Usage

``` r
decomp2sigma(d, G, scale, shape, orientation, ...)
```

## Arguments

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

- ...:

  Catches unused arguments from an indirect or list call via `do.call`.

## Value

A 3-D array whose `[,,k]`th component is the covariance matrix of the
*k*th component in an MVN mixture model.

## See also

[`sigma2decomp`](https://mclust-org.github.io/mclust/reference/sigma2decomp.md)

## Examples

``` r
meEst <- meVEV(iris[,-5], unmap(iris[,5])) 
names(meEst)
#> [1] "modelName"  "prior"      "n"          "d"          "G"         
#> [6] "z"          "parameters" "control"    "loglik"    
meEst$parameters$variance
#> $modelName
#> [1] "VEV"
#> 
#> $d
#> [1] 4
#> 
#> $G
#> [1] 3
#> 
#> $sigma
#> , , 1
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.13321144  0.10939010   0.01919074  0.01158562
#> Sepal.Width    0.10939010  0.15495729   0.01209929  0.01001100
#> Petal.Length   0.01919074  0.01209929   0.02827503  0.00581811
#> Petal.Width    0.01158562  0.01001100   0.00581811  0.01069593
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.22571809  0.07613379   0.14689790  0.04335743
#> Sepal.Width    0.07613379  0.08024536   0.07371576  0.03435377
#> Petal.Length   0.14689790  0.07371576   0.16611738  0.04951832
#> Petal.Width    0.04335743  0.03435377   0.04951832  0.03337836
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.42945746  0.10782916   0.33457182  0.06541484
#> Sepal.Width    0.10782916  0.11595505   0.08903504  0.06133241
#> Petal.Length   0.33457182  0.08903504   0.36427088  0.08711121
#> Petal.Width    0.06541484  0.06133241   0.08711121  0.08666116
#> 
#> 
#> $scale
#> [1] 0.03748690 0.05792051 0.11417100
#> 
#> $shape
#> [1] 6.8586562 1.0115190 0.6265407 0.2300581
#> 
#> $orientation
#> , , 1
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length  -0.66907840   0.5978840   -0.4399628 -0.03607712
#> Sepal.Width   -0.73414783  -0.6206734    0.2746075 -0.01955027
#> Petal.Length  -0.09654390   0.4900556    0.8324495 -0.23990129
#> Petal.Width   -0.06356359   0.1309379    0.1950675  0.96992969
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length    0.7058553  -0.6603954   0.24389224 -0.07850316
#> Sepal.Width     0.3292982   0.5666772   0.68428027  0.31968755
#> Petal.Length    0.5956635   0.3576279  -0.68509739  0.21892660
#> Petal.Width     0.1962550   0.3389057   0.05402459 -0.91853594
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   -0.7104631 -0.30095824   -0.5382163  -0.3391010
#> Sepal.Width    -0.2167218  0.71570645   -0.4318489   0.5042841
#> Petal.Length   -0.6483892 -0.07068137    0.6545951   0.3822313
#> Petal.Width    -0.1669287  0.62625279    0.3087637  -0.6961374
#> 
#> 

dec <- meEst$parameters$variance
decomp2sigma(d=dec$d, G=dec$G, shape=dec$shape, scale=dec$scale,
             orientation = dec$orientation)
#> , , 1
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.13321144 0.10939010 0.01919074 0.01158562
#> [2,] 0.10939010 0.15495729 0.01209929 0.01001100
#> [3,] 0.01919074 0.01209929 0.02827503 0.00581811
#> [4,] 0.01158562 0.01001100 0.00581811 0.01069593
#> 
#> , , 2
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.22571809 0.07613379 0.14689790 0.04335743
#> [2,] 0.07613379 0.08024536 0.07371576 0.03435377
#> [3,] 0.14689790 0.07371576 0.16611738 0.04951832
#> [4,] 0.04335743 0.03435377 0.04951832 0.03337836
#> 
#> , , 3
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.42945746 0.10782916 0.33457182 0.06541484
#> [2,] 0.10782916 0.11595505 0.08903504 0.06133241
#> [3,] 0.33457182 0.08903504 0.36427088 0.08711121
#> [4,] 0.06541484 0.06133241 0.08711121 0.08666116
#> 
#> attr(,"modelName")
#> [1] "VEV"
# \donttest{
do.call("decomp2sigma", dec)  ## alternative call
#> , , 1
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.13321144 0.10939010 0.01919074 0.01158562
#> [2,] 0.10939010 0.15495729 0.01209929 0.01001100
#> [3,] 0.01919074 0.01209929 0.02827503 0.00581811
#> [4,] 0.01158562 0.01001100 0.00581811 0.01069593
#> 
#> , , 2
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.22571809 0.07613379 0.14689790 0.04335743
#> [2,] 0.07613379 0.08024536 0.07371576 0.03435377
#> [3,] 0.14689790 0.07371576 0.16611738 0.04951832
#> [4,] 0.04335743 0.03435377 0.04951832 0.03337836
#> 
#> , , 3
#> 
#>            [,1]       [,2]       [,3]       [,4]
#> [1,] 0.42945746 0.10782916 0.33457182 0.06541484
#> [2,] 0.10782916 0.11595505 0.08903504 0.06133241
#> [3,] 0.33457182 0.08903504 0.36427088 0.08711121
#> [4,] 0.06541484 0.06133241 0.08711121 0.08666116
#> 
#> attr(,"modelName")
#> [1] "VEV"
# }
```
