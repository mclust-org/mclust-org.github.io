# Number of Variance Parameters in Gaussian Mixture Models

Gives the number of variance parameters for parameterizations of the
Gaussian mixture model that are used in MCLUST.

## Usage

``` r
nVarParams(modelName, d, G, ...)
```

## Arguments

- modelName:

  A character string indicating the model. The help file for
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)
  describes the available models.

- d:

  The dimension of the data. Not used for models in which neither the
  shape nor the orientation varies.

- G:

  The number of components in the Gaussian mixture model used to compute
  `loglik`.

- ... :

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

The number of variance parameters in the corresponding Gaussian mixture
model.

## Details

To get the total number of parameters in model, add `G*d` for the means
and `G-1` for the mixing proportions if they are unequal.

## See also

[`bic`](https://mclust-org.github.io/mclust/reference/bic.md),
[`nMclustParams`](https://mclust-org.github.io/mclust/reference/nMclustParams.md).

## Examples

``` r
mapply(nVarParams, mclust.options("emModelNames"), d = 2, G = 3)
#> EII VII EEI VEI EVI VVI EEE VEE EVE VVE EEV VEV EVV VVV 
#>   1   3   2   4   4   6   3   5   5   7   5   7   7   9 
```
