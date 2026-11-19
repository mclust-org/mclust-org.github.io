# Number of Estimated Parameters in Gaussian Mixture Models

Gives the number of estimated parameters for parameterizations of the
Gaussian mixture model that are used in MCLUST.

## Usage

``` r
nMclustParams(modelName, d, G, noise = FALSE, equalPro = FALSE, ...)
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

- noise:

  A logical variable indicating whether or not the model includes an
  optional Poisson noise component.

- equalPro:

  A logical variable indicating whether or not the components in the
  model are assumed to be present in equal proportion.

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
[`nVarParams`](https://mclust-org.github.io/mclust/reference/nVarParams.md).

## Examples

``` r
mapply(nMclustParams, mclust.options("emModelNames"), d = 2, G = 3)
#> EII VII EEI VEI EVI VVI EEE VEE EVE VVE EEV VEV EVV VVV 
#>   9  11  10  12  12  14  11  13  13  15  13  15  15  17 
```
