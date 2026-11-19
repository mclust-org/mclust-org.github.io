# Log-Likelihood of a `MclustDA` object

Returns the log-likelihood for a `MclustDA` object.

## Usage

``` r
# S3 method for class 'MclustDA'
logLik(object, data, ...)
```

## Arguments

- object:

  an object of class `'MclustDA'` resulting from a call to
  [`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

- data:

  the data for which the log-likelihood must be computed. If missing,
  the observed data from the `'MclustDA'` object is used.

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

[`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

## Examples

``` r
# \donttest{
irisMclustDA <- MclustDA(iris[,1:4], iris$Species)
summary(irisMclustDA)
#> ------------------------------------------------ 
#> Gaussian finite mixture model for classification 
#> ------------------------------------------------ 
#> 
#> MclustDA 
#> ------------------------------------------------ 
#>  log-likelihood   n df       BIC
#>       -172.8135 150 49 -591.1482
#>             
#> Classes       n     % Model G
#>   setosa     50 33.33   EEE 2
#>   versicolor 50 33.33   XXX 1
#>   virginica  50 33.33   XXX 1
#> 
#> Training confusion matrix
#> ------------------------------------------------ 
#>             Predicted
#> Classes      setosa versicolor virginica
#>   setosa         50          0         0
#>   versicolor      0         48         2
#>   virginica       0          1        49
#> 
#> Classification error = 0.0200 
#> Brier score          = 0.0116 
logLik(irisMclustDA)
#> 'log Lik.' -172.8135 (df=49)
# }
```
