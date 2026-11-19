# Summary Function for Bootstrap Inference for Gaussian Finite Mixture Models

Summary of bootstrap distribution for the parameters of a Gaussian
mixture model providing either standard errors or percentile bootstrap
confidence intervals.

## Usage

``` r
# S3 method for class 'MclustBootstrap'
summary(object, what = c("se", "ci", "ave"), conf.level = 0.95, ...)
```

## Arguments

- object:

  An object of class `'MclustBootstrap'` as returned by
  [`MclustBootstrap`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md).

- what:

  A character string: `"se"` for the standard errors; `"ci"` for the
  confidence intervals; `"ave"` for the averages.

- conf.level:

  A value specifying the confidence level of the interval.

- ...:

  Further arguments passed to or from other methods.

## Details

For details about the procedure used to obtain the bootstrap
distribution see
[`MclustBootstrap`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md).

## See also

[`MclustBootstrap`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md).

## Examples

``` r
# \donttest{
data(diabetes)
X = diabetes[,-1]
modClust = Mclust(X) 
bootClust = MclustBootstrap(modClust)
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
#> 0.05162109 0.05228003 0.03579922 
#> 
#> Means:
#>                1        2        3
#> glucose 1.084425  3.33190 16.69377
#> insulin 7.876585 27.98349 65.67435
#> sspg    7.876924 31.76700 10.32837
#> 
#> Variances:
#> [,,1]
#>          glucose   insulin      sspg
#> glucose 11.03080  52.39737  51.00843
#> insulin 52.39737 518.32169 413.57951
#> sspg    51.00843 413.57951 640.70793
#> [,,2]
#>           glucose   insulin      sspg
#> glucose  62.97197  627.2149  430.1712
#> insulin 627.21493 7612.4803 3120.8634
#> sspg    430.17120 3120.8634 6872.2195
#> [,,3]
#>           glucose   insulin      sspg
#> glucose  955.0881  3966.236  630.1438
#> insulin 3966.2355 18500.597 2402.5173
#> sspg     630.1438  2402.517  486.3382
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
#> 2.5%  0.4438454 0.1571951 0.1298801
#> 97.5% 0.6500454 0.3581143 0.2683289
#> 
#> Means:
#> [,,1]
#>        glucose  insulin     sspg
#> 2.5%  88.95320 344.9304 150.9777
#> 97.5% 93.14719 375.1357 182.7581
#> [,,2]
#>        glucose  insulin     sspg
#> 2.5%   99.4339 449.4375 263.7769
#> 97.5% 112.8009 561.4841 388.0333
#> [,,3]
#>        glucose   insulin     sspg
#> 2.5%  197.8521  967.8664  60.5601
#> 97.5% 262.1883 1224.9748 101.4817
#> 
#> Variances:
#> [,,1]
#>        glucose  insulin     sspg
#> 2.5%  36.25487 1296.324 1507.570
#> 97.5% 81.50943 3312.735 4172.067
#> [,,2]
#>        glucose   insulin     sspg
#> 2.5%   84.9753  3390.178 12305.56
#> 97.5% 333.2341 32297.036 39317.20
#> [,,3]
#>        glucose   insulin     sspg
#> 2.5%  3453.364  48532.52 1313.878
#> 97.5% 7222.222 122538.59 3184.822

data(acidity)
modDens = densityMclust(acidity, plot = FALSE)
modDens = MclustBootstrap(modDens)
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
#>        1        2 
#> 0.040875 0.040875 
#> 
#> Means:
#>          1          2 
#> 0.04499781 0.06566511 
#> 
#> Variances:
#>          1          2 
#> 0.02362773 0.02362773 
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
#> 2.5%  0.5418831 0.2915891
#> 97.5% 0.7084109 0.4581169
#> 
#> Means:
#>              1        2
#> 2.5%  4.282517 6.184802
#> 97.5% 4.455290 6.448229
#> 
#> Variances:
#>               1         2
#> 2.5%  0.1425153 0.1425153
#> 97.5% 0.2326400 0.2326400
# }
```
