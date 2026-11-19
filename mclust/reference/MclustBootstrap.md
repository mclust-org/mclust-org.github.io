# Resampling-based Inference for Gaussian finite mixture models

Bootstrap or jackknife estimation of standard errors and percentile
bootstrap confidence intervals for the parameters of a Gaussian mixture
model.

## Usage

``` r
MclustBootstrap(object, nboot = 999, 
                type = c("bs", "wlbs", "pb", "jk"),
                alpha = 1, max.nonfit = 10*nboot, 
                verbose = interactive(), ...)
```

## Arguments

- object:

  An object of class `'Mclust'` or `'densityMclust'` providing an
  estimated Gaussian mixture model.

- nboot:

  The number of bootstrap replications.

- type:

  A character string specifying the type of resampling to use:

  `"bs"`

  : nonparametric bootstrap

  `"wlbs"`

  : weighted likelihood bootstrap

  `"pb"`

  : parametric bootstrap

  `"jk"`

  : jackknife

- alpha:

  A numerical value used when `type = "wlbs"` to generate weights from a
  Dirichlet(alpha, ..., alpha) distribution. By default `alpha = 1`, so
  weights are generated from a uniform distribution on the simplex.

- max.nonfit:

  The maximum number of non-estimable models allowed.

- verbose:

  A logical controlling if a text progress bar is displayed during the
  resampling procedure. By default is `TRUE` if the session is
  interactive, and `FALSE` otherwise.

- ...:

  Further arguments passed to or from other methods.

## Details

For a fitted Gaussian mixture model with `object$G` mixture components
and covariances parameterisation `object$modelName`, this function
returns either the bootstrap distribution or the jackknife distribution
of mixture parameters. In the former case, the nonparametric bootstrap
or the weighted likelihood bootstrap approach could be used, so the the
bootstrap procedure generates `nboot` bootstrap samples of the same size
as the original data by resampling with replacement from the observed
data. In the jackknife case, the procedure considers all the samples
obtained by omitting one observation at time.

The resulting resampling distribution can then be used to obtain
standard errors and percentile confidence intervals by the use of
[`summary.MclustBootstrap`](https://mclust-org.github.io/mclust/reference/summary.MclustBootstrap.md)
function.

## Value

An object of class `'MclustBootstrap'` with the following components:

- n:

  The number of observations in the data.

- d:

  The dimension of the data.

- G:

  A value specifying the number of mixture components.

- modelName:

  A character string specifying the mixture model covariances
  parameterisation (see
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)).

- parameters:

  A list of estimated parameters for the mixture components with the
  following components:

  `pro`

  : a vector of mixing proportions.

  `mean`

  : a matrix of means for each component.

  `variance`

  : an array of covariance matrices for each component.

- nboot:

  The number of bootstrap replications if `type = "bs"` or
  `type = "wlbs"`. The sample size if `type = "jk"`.

- type:

  The type of resampling approach used.

- nonfit:

  The number of resamples that did not convergence during the procedure.

- pro:

  A matrix of dimension (`nboot` x `G`) containing the bootstrap
  distribution for the mixing proportion.

- mean:

  An array of dimension (`nboot` x `d` x `G`), where `d` is the
  dimension of the data, containing the bootstrap distribution for the
  component means.

- variance:

  An array of dimension (`nboot` x `d` x `d` x `G`), where `d` is the
  dimension of the data, containing the bootstrap distribution for the
  component covariances.

## References

Davison, A. and Hinkley, D. (1997) *Bootstrap Methods and Their
Applications*. Cambridge University Press.

McLachlan, G.J. and Peel, D. (2000) *Finite Mixture Models*. Wiley.

O'Hagan A., Murphy T. B., Gormley I. C. and Scrucca L. (2015) On
Estimation of Parameter Uncertainty in Model-Based Clustering. Submitted
to *Computational Statistics*.

## See also

[`summary.MclustBootstrap`](https://mclust-org.github.io/mclust/reference/summary.MclustBootstrap.md),
[`plot.MclustBootstrap`](https://mclust-org.github.io/mclust/reference/plot.MclustBoostrap.md),
[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md),
[`densityMclust`](https://mclust-org.github.io/mclust/reference/densityMclust.md).

## Examples

``` r
# \donttest{
data(diabetes)
X <- diabetes[,-1]
modClust <- Mclust(X) 
bootClust <- MclustBootstrap(modClust)
summary(bootClust, what = "se")
#> ---------------------------------------------------------- 
#> Resampling standard errors 
#> ---------------------------------------------------------- 
#> Model                      = VVV 
#> Num. of mixture components = 3 
#> Replications               = 999 
#> Type                       = nonparametric bootstrap 
#> 
#> Mixing probabilities:
#>          1          2          3 
#> 0.05233913 0.04934156 0.03494921 
#> 
#> Means:
#>                1         2         3
#> glucose 1.006765  3.292138 16.647842
#> insulin 7.843462 27.548758 65.966648
#> sspg    7.779137 30.692950  9.983326
#> 
#> Variances:
#> [,,1]
#>          glucose   insulin      sspg
#> glucose 10.64753  50.63956  48.61746
#> insulin 50.63956 524.57448 404.58548
#> sspg    48.61746 404.58548 647.28663
#> [,,2]
#>           glucose   insulin      sspg
#> glucose  63.32924  615.6663  435.9384
#> insulin 615.66626 7336.6348 3086.8630
#> sspg    435.93842 3086.8630 6780.6833
#> [,,3]
#>           glucose   insulin      sspg
#> glucose  995.5757  4072.436  636.6375
#> insulin 4072.4359 18497.223 2422.3585
#> sspg     636.6375  2422.358  473.1342
summary(bootClust, what = "ci")
#> ---------------------------------------------------------- 
#> Resampling confidence intervals 
#> ---------------------------------------------------------- 
#> Model                      = VVV 
#> Num. of mixture components = 3 
#> Replications               = 999 
#> Type                       = nonparametric bootstrap 
#> Confidence level           = 0.95 
#> 
#> Mixing probabilities:
#>               1         2         3
#> 2.5%  0.4478146 0.1543735 0.1366100
#> 97.5% 0.6435917 0.3555992 0.2706844
#> 
#> Means:
#> [,,1]
#>        glucose  insulin     sspg
#> 2.5%  89.07430 344.1764 150.6252
#> 97.5% 93.22599 375.1545 182.3419
#> [,,2]
#>         glucose  insulin     sspg
#> 2.5%   98.94194 452.8820 258.8573
#> 97.5% 112.41549 558.5595 382.9069
#> [,,3]
#>        glucose   insulin      sspg
#> 2.5%  195.9121  962.2439  61.13393
#> 97.5% 263.0002 1228.3581 101.71555
#> 
#> Variances:
#> [,,1]
#>        glucose  insulin     sspg
#> 2.5%  38.64074 1257.854 1493.575
#> 97.5% 80.71107 3239.885 4053.007
#> [,,2]
#>         glucose   insulin     sspg
#> 2.5%   91.03058  3405.421 12808.23
#> 97.5% 345.01024 30217.896 38356.71
#> [,,3]
#>        glucose   insulin     sspg
#> 2.5%  3437.032  45919.79 1392.434
#> 97.5% 7452.134 119570.85 3311.401

data(acidity)
modDens <- densityMclust(acidity, plot = FALSE)
modDens <- MclustBootstrap(modDens)
summary(modDens, what = "se")
#> ---------------------------------------------------------- 
#> Resampling standard errors 
#> ---------------------------------------------------------- 
#> Model                      = E 
#> Num. of mixture components = 2 
#> Replications               = 999 
#> Type                       = nonparametric bootstrap 
#> 
#> Mixing probabilities:
#>          1          2 
#> 0.04171062 0.04171062 
#> 
#> Means:
#>          1          2 
#> 0.04515616 0.06660559 
#> 
#> Variances:
#>          1          2 
#> 0.02362844 0.02362844 
summary(modDens, what = "ci")
#> ---------------------------------------------------------- 
#> Resampling confidence intervals 
#> ---------------------------------------------------------- 
#> Model                      = E 
#> Num. of mixture components = 2 
#> Replications               = 999 
#> Type                       = nonparametric bootstrap 
#> Confidence level           = 0.95 
#> 
#> Mixing probabilities:
#>               1         2
#> 2.5%  0.5388678 0.2989128
#> 97.5% 0.7010872 0.4611322
#> 
#> Means:
#>              1        2
#> 2.5%  4.285496 6.179751
#> 97.5% 4.458050 6.442143
#> 
#> Variances:
#>               1         2
#> 2.5%  0.1390296 0.1390296
#> 97.5% 0.2297464 0.2297464
# }
```
