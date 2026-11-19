# EM algorithm with weights starting with M-step for parameterized Gaussian mixture models

Implements the EM algorithm for fitting Gaussian mixture models
parameterized by eigenvalue decomposition, when observations have
weights, starting with the maximization step.

## Usage

``` r
me.weighted(data, modelName, z, weights = NULL, prior = NULL, 
            control = emControl(), Vinv = NULL, warn = NULL, ...)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- modelName:

  A character string indicating the model. The help file for
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)
  describes the available models.

- z:

  A matrix whose `[i,k]`th entry is an initial estimate of the
  conditional probability of the ith observation belonging to the *k*th
  component of the mixture.

- weights:

  A vector of positive weights, where the `[i]`th entry is the weight
  for the ith observation. If any of the weights are greater than one,
  then they are scaled so that the maximum weight is one.

- prior:

  Specification of a conjugate prior on the means and variances. See the
  help file for `priorControl` for further information. The default
  assumes no prior.

- control:

  A list of control parameters for EM. The defaults are set by the call
  [`emControl`](https://mclust-org.github.io/mclust/reference/emControl.md).

- Vinv:

  If the model is to include a noise term, `Vinv` is an estimate of the
  reciprocal hypervolume of the data region. If set to a negative value
  or 0, the model will include a noise term with the reciprocal
  hypervolume estimated by the function `hypvol`. The default is not to
  assume a noise term in the model through the setting `Vinv=NULL`.

- warn:

  A logical value indicating whether or not certain warnings (usually
  related to singularity) should be issued when the estimation fails.
  The default is set by `warn` using
  [`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md).

- ...:

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A list including the following components:

- modelName:

  A character string identifying the model (same as the input argument).

- z:

  A matrix whose `[i,k]`th entry is the conditional probability of the
  *i*th observation belonging to the *k*th component of the mixture.

- parameters:

  `pro`

  : A vector whose *k*th component is the mixing proportion for the
    *k*th component of the mixture model. If the model includes a
    Poisson term for noise, there should be one more mixing proportion
    than the number of Gaussian components.

  `mean`

  : The mean for each component. If there is more than one component,
    this is a matrix whose kth column is the mean of the *k*th component
    of the mixture model.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

  `Vinv`

  : The estimate of the reciprocal hypervolume of the data region used
    in the computation when the input indicates the addition of a noise
    component to the model.

- loglik:

  The log-likelihood for the estimated mixture model.

- bic:

  The BIC value for the estimated mixture model.

- Attributes::

  `"info"` Information on the iteration.\
  `"WARNING"` An appropriate warning if problems are encountered in the
  computations.

## Details

This is a more efficient version made available with mclust \\(ge 6.1\\)
using Fortran code internally.

## Author

T. Brendan Murphy, Luca Scrucca

## See also

[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`meE`](https://mclust-org.github.io/mclust/reference/meE.md), ...,
[`meVVV`](https://mclust-org.github.io/mclust/reference/meE.md),
[`em`](https://mclust-org.github.io/mclust/reference/em.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`priorControl`](https://mclust-org.github.io/mclust/reference/priorControl.md),
[`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md),
[`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md)

## Examples

``` r
w = rexp(nrow(iris))
w = w/mean(w)
c(summary(w), sum = sum(w))
#>         Min.      1st Qu.       Median         Mean      3rd Qu.         Max. 
#>   0.00206268   0.32520507   0.71293241   1.00000000   1.47179878   6.53900532 
#>          sum 
#> 150.00000000 
z = unmap(sample(1:3, size = nrow(iris), replace = TRUE))
MEW = me.weighted(data = iris[,-5], modelName = "VVV", 
                  z = z, weights = w)
str(MEW,1)
#> List of 10
#>  $ modelName : chr "VVV"
#>  $ prior     : NULL
#>  $ n         : int 150
#>  $ d         : int 4
#>  $ G         : int 3
#>  $ z         : num [1:150, 1:3] 2.89e-18 6.45e-12 7.86e-15 3.23e-12 2.31e-19 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ parameters:List of 3
#>  $ weights   : num [1:150] 0.0267 0.0986 0.2817 0.029 0.429 ...
#>  $ loglik    : num -178
#>  $ bic       : num -576
#>  - attr(*, "returnCode")= num 0
```
