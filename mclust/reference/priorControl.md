# Conjugate Prior for Gaussian Mixtures.

Specify a conjugate prior for Gaussian mixtures.

## Usage

``` r
priorControl(functionName = "defaultPrior", ...)
```

## Arguments

- functionName:

  The name of the function specifying the conjugate prior. By default
  the function
  [`defaultPrior`](https://mclust-org.github.io/mclust/reference/defaultPrior.md)
  is used, and this can also be used as a template for alternative
  specification.

- ...:

  Optional named arguments to the function specified in `functionName`
  together with their values.

## Value

A list with the function name as the first component. The remaining
components (if any) consist of a list of arguments to the function with
assigned values.

## Details

The function `priorControl` is used to specify a conjugate prior for EM
within *MCLUST*.\
Note that, as described in
[`defaultPrior`](https://mclust-org.github.io/mclust/reference/defaultPrior.md),
in the multivariate case only 10 out of 14 models may be used in
conjunction with a prior, i.e. those available in *MCLUST* up to version
4.4.

## References

C. Fraley and A. E. Raftery (2007). Bayesian regularization for normal
mixture estimation and model-based clustering. *Journal of
Classification 24:155-181*.

## See also

[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`defaultPrior`](https://mclust-org.github.io/mclust/reference/defaultPrior.md)

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
```
