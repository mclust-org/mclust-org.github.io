# Component Density for Parameterized MVN Mixture Models

Computes component densities for observations in MVN mixture models
parameterized by eigenvalue decomposition.

## Usage

``` r
cdens(data, modelName, parameters, logarithm = FALSE, warn = NULL, ...)
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

- parameters:

  The parameters of the model:

  `mean`

  : The mean for each component. If there is more than one component,
    this is a matrix whose kth column is the mean of the *k*th component
    of the mixture model.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

- logarithm:

  A logical value indicating whether or not the logarithm of the
  component densities should be returned. The default is to return the
  component densities, obtained from the log component densities by
  exponentiation.

- warn:

  A logical value indicating whether or not a warning should be issued
  when computations fail. The default is `warn=FALSE`.

- ...:

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A numeric matrix whose `[i,k]`th entry is the density or log density of
observation *i* in component *k*. The densities are not scaled by mixing
proportions.

## Note

When one or more component densities are very large in magnitude, it may
be possible to compute the logarithm of the component densities but not
the component densities themselves due to overflow.

## See also

[`cdensE`](https://mclust-org.github.io/mclust/reference/cdensE.md),
...,
[`cdensVVV`](https://mclust-org.github.io/mclust/reference/cdensE.md),
[`dens`](https://mclust-org.github.io/mclust/reference/dens.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md),
[`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md),
[`do.call`](https://rdrr.io/r/base/do.call.html)

## Examples

``` r
z2 <- unmap(hclass(hcVVV(faithful),2)) # initial value for 2 class case

model <- me(modelName = "EEE", data = faithful, z = z2)
cdens(modelName = "EEE", data = faithful, logarithm = TRUE, 
      parameters = model$parameters)[1:5,]
#>         [,1]       [,2]
#> 1  -4.504988 -15.570952
#> 2 -28.297451  -2.782971
#> 3  -6.042371 -11.151408
#> 4 -18.516948  -3.349049
#> 5  -2.967399 -30.134605

data(cross)
odd <- seq(1, nrow(cross), by = 2)
oddBIC <- mclustBIC(cross[odd,-1]) 
oddModel <- mclustModel(cross[odd,-1], oddBIC) ## best parameter estimates
names(oddModel)
#> [1] "modelName"  "n"          "d"          "G"          "bic"       
#> [6] "loglik"     "parameters" "z"         

even <- odd + 1
densities <- cdens(modelName = oddModel$modelName, data = cross[even,-1], 
                   parameters = oddModel$parameters)
cbind(class = cross[even,1], densities)[1:5,]
#>    class                         
#> 2      1 0.014497007 1.123598e-07
#> 4      1 0.001087501 1.333448e-75
#> 6      1 0.005431090 1.180715e-02
#> 8      1 0.013545695 2.488475e-10
#> 10     1 0.001215779 9.250963e-03
```
