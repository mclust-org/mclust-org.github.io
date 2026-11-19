# Univariate or Multivariate Normal Fit

Computes the mean, covariance, and log-likelihood from fitting a single
Gaussian to given data (univariate or multivariate normal).

## Usage

``` r
mvn( modelName, data, prior = NULL, warn = NULL, ...)
```

## Arguments

- modelName:

  A character string representing a model name. This can be either
  `"Spherical"`, `"Diagonal"`, or `"Ellipsoidal"` or else\
  `"X"` for one-dimensional data,\
  `"XII"` for a spherical Gaussian,\
  `"XXI"` for a diagonal Gaussian\
  `"XXX"` for a general ellipsoidal Gaussian

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- prior:

  Specification of a conjugate prior on the means and variances. The
  default assumes no prior.

- warn:

  A logical value indicating whether or not a warning should be issued
  whenever a singularity is encountered. The default is given by
  `mclust.options("warn")`.

- ... :

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A list including the following components:

- modelName:

  A character string identifying the model (same as the input argument).

- parameters:

  `mean`

  : The mean for each component. If there is more than one component,
    this is a matrix whose kth column is the mean of the *k*th component
    of the mixture model.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

- loglik:

  The log likelihood for the data in the mixture model.

- Attributes::

  `"WARNING"` An appropriate warning if problems are encountered in the
  computations.

## See also

[`mvnX`](https://mclust-org.github.io/mclust/reference/mvnX.md),
[`mvnXII`](https://mclust-org.github.io/mclust/reference/mvnX.md),
[`mvnXXI`](https://mclust-org.github.io/mclust/reference/mvnX.md),
[`mvnXXX`](https://mclust-org.github.io/mclust/reference/mvnX.md),
[`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)

## Examples

``` r
n <- 1000

set.seed(0)
x <- rnorm(n, mean = -1, sd = 2)
mvn(modelName = "X", x) 
#> $modelName
#> [1] "X"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 1
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#> [1] -1.031659
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "X"
#> 
#> $parameters$variance$d
#> [1] 1
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$sigmasq
#> [1] 3.980071
#> 
#> 
#> 
#> $loglik
#> [1] -2109.588
#> 
#> attr(,"returnCode")
#> [1] 0

mu <- c(-1, 0, 1)

set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) %*% (2*diag(3)), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvn(modelName = "XII", x) 
#> $modelName
#> [1] "XII"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 3
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#>             [,1]
#> [1,] -1.03165915
#> [2,] -0.04957289
#> [3,]  1.13628546
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "XII"
#> 
#> $parameters$variance$d
#> [1] 3
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$sigmasq
#> [1] 3.963291
#> 
#> $parameters$variance$Sigma
#>          [,1]     [,2]     [,3]
#> [1,] 3.963291 0.000000 0.000000
#> [2,] 0.000000 3.963291 0.000000
#> [3,] 0.000000 0.000000 3.963291
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>          [,1]     [,2]     [,3]
#> [1,] 3.963291 0.000000 0.000000
#> [2,] 0.000000 3.963291 0.000000
#> [3,] 0.000000 0.000000 3.963291
#> 
#> 
#> $parameters$variance$scale
#> [1] 3.963291
#> 
#> 
#> 
#> $loglik
#> [1] -6322.428
#> 
#> attr(,"returnCode")
#> [1] 0
mvn(modelName = "Spherical", x) 
#> $modelName
#> [1] "XII"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 3
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#>             [,1]
#> [1,] -1.03165915
#> [2,] -0.04957289
#> [3,]  1.13628546
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "XII"
#> 
#> $parameters$variance$d
#> [1] 3
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$sigmasq
#> [1] 3.963291
#> 
#> $parameters$variance$Sigma
#>          [,1]     [,2]     [,3]
#> [1,] 3.963291 0.000000 0.000000
#> [2,] 0.000000 3.963291 0.000000
#> [3,] 0.000000 0.000000 3.963291
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>          [,1]     [,2]     [,3]
#> [1,] 3.963291 0.000000 0.000000
#> [2,] 0.000000 3.963291 0.000000
#> [3,] 0.000000 0.000000 3.963291
#> 
#> 
#> $parameters$variance$scale
#> [1] 3.963291
#> 
#> 
#> 
#> $loglik
#> [1] -6322.428
#> 
#> attr(,"returnCode")
#> [1] 0

set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) %*% diag(1:3), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvn(modelName = "XXI", x)
#> $modelName
#> [1] "XXI"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 3
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#>             [,1]
#> [1,] -1.01582957
#> [2,] -0.04957289
#> [3,]  1.20442820
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "XXI"
#> 
#> $parameters$variance$d
#> [1] 3
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$Sigma
#>           [,1]    [,2]     [,3]
#> [1,] 0.9950176 0.00000 0.000000
#> [2,] 0.0000000 4.27032 0.000000
#> [3,] 0.0000000 0.00000 8.188836
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>           [,1]    [,2]     [,3]
#> [1,] 0.9950176 0.00000 0.000000
#> [2,] 0.0000000 4.27032 0.000000
#> [3,] 0.0000000 0.00000 8.188836
#> 
#> 
#> $parameters$variance$scale
#> [1] 3.264659
#> 
#> $parameters$variance$shape
#> [1] 0.3047846 1.3080448 2.5083282
#> 
#> 
#> 
#> $loglik
#> [1] -6031.548
#> 
#> attr(,"returnCode")
#> [1] 0
mvn(modelName = "Diagonal", x)
#> $modelName
#> [1] "XXI"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 3
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#>             [,1]
#> [1,] -1.01582957
#> [2,] -0.04957289
#> [3,]  1.20442820
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "XXI"
#> 
#> $parameters$variance$d
#> [1] 3
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$Sigma
#>           [,1]    [,2]     [,3]
#> [1,] 0.9950176 0.00000 0.000000
#> [2,] 0.0000000 4.27032 0.000000
#> [3,] 0.0000000 0.00000 8.188836
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>           [,1]    [,2]     [,3]
#> [1,] 0.9950176 0.00000 0.000000
#> [2,] 0.0000000 4.27032 0.000000
#> [3,] 0.0000000 0.00000 8.188836
#> 
#> 
#> $parameters$variance$scale
#> [1] 3.264659
#> 
#> $parameters$variance$shape
#> [1] 0.3047846 1.3080448 2.5083282
#> 
#> 
#> 
#> $loglik
#> [1] -6031.548
#> 
#> attr(,"returnCode")
#> [1] 0

Sigma <- matrix(c(9,-4,1,-4,9,4,1,4,9), 3, 3)
set.seed(0)
x <- sweep(matrix(rnorm(n*3), n, 3) %*% chol(Sigma), 
           MARGIN = 2, STATS = mu, FUN = "+")
mvn(modelName = "XXX", x) 
#> $modelName
#> [1] "XXX"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 3
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#>             [,1]
#> [1,] -1.04748872
#> [2,] -0.04550547
#> [3,]  1.12277306
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "XXX"
#> 
#> $parameters$variance$d
#> [1] 3
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$Sigma
#>           [,1]      [,2]     [,3]
#> [1,]  8.955159 -4.084953 1.015574
#> [2,] -4.084953  9.572449 4.417112
#> [3,]  1.015574  4.417112 8.800769
#> 
#> $parameters$variance$cholSigma
#>          [,1]      [,2]       [,3]
#> [1,] 2.992517 -1.365056  0.3393713
#> [2,] 0.000000 -2.776521 -1.7577291
#> [3,] 0.000000  0.000000 -2.3655833
#> 
#> $parameters$variance$cholsigma
#>          [,1]      [,2]       [,3]
#> [1,] 2.992517 -1.365056  0.3393713
#> [2,] 0.000000 -2.776521 -1.7577291
#> [3,] 0.000000  0.000000 -2.3655833
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>           [,1]      [,2]     [,3]
#> [1,]  8.955159 -4.084953 1.015574
#> [2,] -4.084953  9.572449 4.417112
#> [3,]  1.015574  4.417112 8.800769
#> 
#> 
#> 
#> 
#> $loglik
#> [1] -7235.154
#> 
#> attr(,"returnCode")
#> [1] 0
mvn(modelName = "Ellipsoidal", x) 
#> $modelName
#> [1] "XXX"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 1000
#> 
#> $d
#> [1] 3
#> 
#> $G
#> [1] 1
#> 
#> $parameters
#> $parameters$pro
#> [1] 1
#> 
#> $parameters$mean
#>             [,1]
#> [1,] -1.04748872
#> [2,] -0.04550547
#> [3,]  1.12277306
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "XXX"
#> 
#> $parameters$variance$d
#> [1] 3
#> 
#> $parameters$variance$G
#> [1] 1
#> 
#> $parameters$variance$Sigma
#>           [,1]      [,2]     [,3]
#> [1,]  8.955159 -4.084953 1.015574
#> [2,] -4.084953  9.572449 4.417112
#> [3,]  1.015574  4.417112 8.800769
#> 
#> $parameters$variance$cholSigma
#>          [,1]      [,2]       [,3]
#> [1,] 2.992517 -1.365056  0.3393713
#> [2,] 0.000000 -2.776521 -1.7577291
#> [3,] 0.000000  0.000000 -2.3655833
#> 
#> $parameters$variance$cholsigma
#>          [,1]      [,2]       [,3]
#> [1,] 2.992517 -1.365056  0.3393713
#> [2,] 0.000000 -2.776521 -1.7577291
#> [3,] 0.000000  0.000000 -2.3655833
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>           [,1]      [,2]     [,3]
#> [1,]  8.955159 -4.084953 1.015574
#> [2,] -4.084953  9.572449 4.417112
#> [3,]  1.015574  4.417112 8.800769
#> 
#> 
#> 
#> 
#> $loglik
#> [1] -7235.154
#> 
#> attr(,"returnCode")
#> [1] 0
```
