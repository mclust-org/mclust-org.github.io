# Default conjugate prior for Gaussian mixtures

Default conjugate prior specification for Gaussian mixtures.

## Usage

``` r
defaultPrior(data, G, modelName, ...)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- G:

  The number of mixture components.

- modelName:

  A character string indicating the model:\
  \
  `"E"`: equal variance (univariate)\
  `"V"`: variable variance (univariate)\
  `"EII"`: spherical, equal volume\
  `"VII"`: spherical, unequal volume\
  `"EEI"`: diagonal, equal volume and shape\
  `"VEI"`: diagonal, varying volume, equal shape\
  `"EVI"`: diagonal, equal volume, varying shape\
  `"VVI"`: diagonal, varying volume and shape\
  `"EEE"`: ellipsoidal, equal volume, shape, and orientation\
  `"EEV"`: ellipsoidal, equal volume and equal shape\
  `"VEV"`: ellipsoidal, equal shape\
  `"VVV"`: ellipsoidal, varying volume, shape, and orientation.\
  \
  A description of the models above is provided in the help of
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md).
  Note that in the multivariate case only 10 out of 14 models may be
  used in conjunction with a prior, i.e. those available in *MCLUST* up
  to version 4.4.

- ...:

  One or more of the following:

  `dof`

  : The degrees of freedom for the prior on the variance. The default is
    `d + 2`, where `d` is the dimension of the data.

  `scale`

  : The scale parameter for the prior on the variance. The default is
    `var(data)/G^(2/d)`, where `d` is the dimension of the data.

  `shrinkage`

  : The shrinkage parameter for the prior on the mean. The default value
    is 0.01. If 0 or NA, no prior is assumed for the mean.

  `mean`

  : The mean parameter for the prior. The default value is
    `colMeans(data)`.

## Value

A list giving the prior degrees of freedom, scale, shrinkage, and mean.

## Details

`defaultPrior` is a function whose default is to output the default
prior specification for EM within *MCLUST*.\
Furthermore, `defaultPrior` can be used as a template to specify
alternative parameters for a conjugate prior.

## References

C. Fraley and A. E. Raftery (2002). Model-based clustering, discriminant
analysis, and density estimation. *Journal of the American Statistical
Association* 97:611-631.

C. Fraley and A. E. Raftery (2005, revised 2009). Bayesian
regularization for normal mixture estimation and model-based clustering.
Technical Report, Department of Statistics, University of Washington.

C. Fraley and A. E. Raftery (2007). Bayesian regularization for normal
mixture estimation and model-based clustering. *Journal of
Classification* 24:155-181.

## See also

[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`priorControl`](https://mclust-org.github.io/mclust/reference/priorControl.md)

## Examples

``` r
# default prior
irisBIC <- mclustBIC(iris[,-5], prior = priorControl())
#> Warning: The presence of BIC values equal to NA is likely due to one or more of the mixture proportions being estimated as zero, so that the model estimated reduces to one with a smaller number of components.
summary(irisBIC, iris[,-5])
#> Best BIC values:
#>              VEV,2       VEV,3      VVV,2
#> BIC      -580.8136 -587.403843 -592.51283
#> BIC diff    0.0000   -6.590289  -11.69928
#> 
#> Classification table for model (VEV,2): 
#> 
#>   1   2 
#>  50 100 

# equivalent to previous example
irisBIC <- mclustBIC(iris[,-5], 
                     prior = priorControl(functionName = "defaultPrior"))
#> Warning: The presence of BIC values equal to NA is likely due to one or more of the mixture proportions being estimated as zero, so that the model estimated reduces to one with a smaller number of components.
summary(irisBIC, iris[,-5])
#> Best BIC values:
#>              VEV,2       VEV,3      VVV,2
#> BIC      -580.8136 -587.403843 -592.51283
#> BIC diff    0.0000   -6.590289  -11.69928
#> 
#> Classification table for model (VEV,2): 
#> 
#>   1   2 
#>  50 100 

# no prior on the mean; default prior on variance
irisBIC <- mclustBIC(iris[,-5], prior = priorControl(shrinkage = 0))
#> Warning: The presence of BIC values equal to NA is likely due to one or more of the mixture proportions being estimated as zero, so that the model estimated reduces to one with a smaller number of components.
summary(irisBIC, iris[,-5])
#> Best BIC values:
#>              VEV,2       VEV,3      VVV,2
#> BIC      -580.2861 -586.792195 -592.07132
#> BIC diff    0.0000   -6.506112  -11.78523
#> 
#> Classification table for model (VEV,2): 
#> 
#>   1   2 
#>  50 100 

# equivalent to previous example
irisBIC <- mclustBIC(iris[,-5], prior =
                     priorControl(functionName="defaultPrior", shrinkage=0))
#> Warning: The presence of BIC values equal to NA is likely due to one or more of the mixture proportions being estimated as zero, so that the model estimated reduces to one with a smaller number of components.
summary(irisBIC, iris[,-5])
#> Best BIC values:
#>              VEV,2       VEV,3      VVV,2
#> BIC      -580.2861 -586.792195 -592.07132
#> BIC diff    0.0000   -6.506112  -11.78523
#> 
#> Classification table for model (VEV,2): 
#> 
#>   1   2 
#>  50 100 

defaultPrior( iris[-5], G = 3, modelName = "VVV")
#> $shrinkage
#> [1] 0.01
#> 
#> $mean
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>     5.843333     3.057333     3.758000     1.199333 
#> 
#> $dof
#> [1] 6
#> 
#> $scale
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.39588533 -0.02449928    0.7357264  0.29806902
#> Sepal.Width   -0.02449928  0.10968467   -0.1903272 -0.07022853
#> Petal.Length   0.73572636 -0.19032720    1.7991839  0.74802043
#> Petal.Width    0.29806902 -0.07022853    0.7480204  0.33544412
#> 
```
