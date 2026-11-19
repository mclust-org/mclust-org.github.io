# Log-Likelihood of a `Mclust` object

Returns the log-likelihood for a `'Mclust'` object.

## Usage

``` r
# S3 method for class 'Mclust'
logLik(object, ...)
```

## Arguments

- object:

  an object of class `'Mclust'` resulting from a call to
  [`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md).

- ...:

  further arguments passed to or from other methods.

## Value

Returns an object of class `'logLik'` with an element providing the
maximized log-likelihood, and further arguments giving the number of
(estimated) parameters in the model (`"df"`) and the sample size
(`"nobs"`).

## Author

Luca Scrucca

## See also

[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md).

## Examples

``` r
# \donttest{
irisMclust <- Mclust(iris[,1:4])
summary(irisMclust)
#> ---------------------------------------------------- 
#> Gaussian finite mixture model fitted by EM algorithm 
#> ---------------------------------------------------- 
#> 
#> Mclust VEV (ellipsoidal, equal shape) model with 2 components: 
#> 
#>  log-likelihood   n df       BIC       ICL
#>        -215.726 150 26 -561.7285 -561.7289
#> 
#> Clustering table:
#>   1   2 
#>  50 100 
logLik(irisMclust)
#> 'log Lik.' -215.726 (df=26)
# }
```
