# BIC for Parameterized Gaussian Mixture Models

Computes the BIC (Bayesian Information Criterion) for parameterized
mixture models given the loglikelihood, the dimension of the data, and
number of mixture components in the model.

## Usage

``` r
bic(modelName, loglik, n, d, G, noise=FALSE, equalPro=FALSE, ...)
```

## Arguments

- modelName:

  A character string indicating the model. The help file for
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)
  describes the available models.

- loglik:

  The log-likelihood for a data set with respect to the Gaussian mixture
  model specified in the `modelName` argument.

- n:

  The number of observations in the data used to compute `loglik`.

- d:

  The dimension of the data used to compute `loglik`.

- G:

  The number of components in the Gaussian mixture model used to compute
  `loglik`.

- noise:

  A logical variable indicating whether or not the model includes an
  optional Poisson noise component. The default is to assume no noise
  component.

- equalPro:

  A logical variable indicating whether or not the components in the
  model are assumed to be present in equal proportion. The default is to
  assume unequal mixing proportions.

- ...:

  Catches unused arguments in an indirect or list call via `do.call`.

## Value

The BIC or Bayesian Information Criterion for the given input arguments.

## See also

[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md),
[`nVarParams`](https://mclust-org.github.io/mclust/reference/nVarParams.md),
[`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md).

## Examples

``` r
# \donttest{
n <- nrow(iris)
d <- ncol(iris)-1
G <- 3

emEst <- me(modelName="VVI", data=iris[,-5], unmap(iris[,5]))
names(emEst)
#> [1] "modelName"  "prior"      "n"          "d"          "G"         
#> [6] "z"          "parameters" "control"    "loglik"    

args(bic)
#> function (modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE, 
#>     ...) 
#> NULL
bic(modelName="VVI", loglik=emEst$loglik, n=n, d=d, G=G)
#> [1] -744.0273
# do.call("bic", emEst)    ## alternative call
# }
```
