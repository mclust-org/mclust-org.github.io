# Summary function for model-based clustering via BIC

Optimal model characteristics and classification for model-based
clustering via `mclustBIC`.

## Usage

``` r
# S3 method for class 'mclustBIC'
summary(object, data, G, modelNames, ...)
```

## Arguments

- object:

  An `'mclustBIC'` object, which is the result of applying `mclustBIC`
  to `data`.

- data:

  The matrix or vector of observations used to generate \`object'.

- G:

  A vector of integers giving the numbers of mixture components
  (clusters) from which the best model according to BIC will be selected
  (`as.character(G)` must be a subset of the row names of `object`). The
  default is to select the best model for all numbers of mixture
  components used to obtain `object`.

- modelNames:

  A vector of integers giving the model parameterizations from which the
  best model according to BIC will be selected (`as.character(model)`
  must be a subset of the column names of `object`). The default is to
  select the best model for parameterizations used to obtain `object`.

- ...:

  Not used. For generic/method consistency.

## Value

A list giving the optimal (according to BIC) parameters, conditional
probabilities `z`, and log-likelihood, together with the associated
classification and its uncertainty.

The details of the output components are as follows:

- modelName:

  A character string denoting the model corresponding to the optimal
  BIC.

- n:

  The number of observations in the data.

- d:

  The dimension of the data.

- G:

  The number of mixture components in the model corresponding to the
  optimal BIC.

- bic:

  The optimal BIC value.

- loglik:

  The log-likelihood corresponding to the optimal BIC.

- parameters:

  A list with the following components:

  `pro`

  : A vector whose *k*th component is the mixing proportion for the
    *k*th component of the mixture model. If missing, equal proportions
    are assumed.

  `mean`

  : The mean for each component. If there is more than one component,
    this is a matrix whose kth column is the mean of the *k*th component
    of the mixture model.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

- z:

  A matrix whose *\[i,k\]*th entry is the probability that observation
  *i* in the data belongs to the *k*th class.

- classification:

  `map(z)`: The classification corresponding to `z`.

- uncertainty:

  The uncertainty associated with the classification.

- Attributes::

  `"bestBICvalues"` Some of the best bic values for the analysis.\
  `"prior"` The prior as specified in the input.\
  `"control"` The control parameters for EM as specified in the input.\
  `"initialization"` The parameters used to initial EM for computing the
  maximum likelihood values used to obtain the BIC.

## See also

[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)
[`mclustModel`](https://mclust-org.github.io/mclust/reference/mclustModel.md)

## Examples

``` r
irisBIC <- mclustBIC(iris[,-5])
summary(irisBIC, iris[,-5])
#> Best BIC values:
#>              VEV,2        VEV,3      VVV,2
#> BIC      -561.7285 -562.5522369 -574.01783
#> BIC diff    0.0000   -0.8237748  -12.28937
#> 
#> Classification table for model (VEV,2): 
#> 
#>   1   2 
#>  50 100 
summary(irisBIC, iris[,-5], G = 1:6, modelNames = c("VII", "VVI", "VVV"))
#> Best BIC values:
#>              VVV,2       VVV,3      VVV,4
#> BIC      -574.0178 -580.839630 -630.59996
#> BIC diff    0.0000   -6.821798  -56.58213
#> 
#> Classification table for model (VVV,2): 
#> 
#>   1   2 
#>  50 100 
```
