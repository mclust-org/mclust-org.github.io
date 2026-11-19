# Summarizing Gaussian Finite Mixture Model Fits

Summary method for class `"Mclust"`.

## Usage

``` r
# S3 method for class 'Mclust'
summary(object, classification = TRUE, parameters = FALSE, ...)
# S3 method for class 'summary.Mclust'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- object:

  An object of class `'Mclust'` resulting of a call to
  [`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md) or
  [`densityMclust`](https://mclust-org.github.io/mclust/reference/densityMclust.md).

- x:

  An object of class `'summary.Mclust'`, usually, a result of a call to
  `summary.Mclust`.

- classification:

  Logical; if `TRUE` a table of MAP classification/clustering of
  observations is printed.

- parameters:

  Logical; if `TRUE`, the parameters of mixture components are printed.

- digits:

  The number of significant digits to use when printing.

- ...:

  Further arguments passed to or from other methods.

## Author

Luca Scrucca

## See also

[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md),
[`densityMclust`](https://mclust-org.github.io/mclust/reference/densityMclust.md).

## Examples

``` r
# \donttest{
mod1 = Mclust(iris[,1:4])
summary(mod1)
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
summary(mod1, parameters = TRUE, classification = FALSE)
#> ---------------------------------------------------- 
#> Gaussian finite mixture model fitted by EM algorithm 
#> ---------------------------------------------------- 
#> 
#> Mclust VEV (ellipsoidal, equal shape) model with 2 components: 
#> 
#>  log-likelihood   n df       BIC       ICL
#>        -215.726 150 26 -561.7285 -561.7289
#> 
#> Mixing probabilities:
#>         1         2 
#> 0.3333319 0.6666681 
#> 
#> Means:
#>                   [,1]     [,2]
#> Sepal.Length 5.0060022 6.261996
#> Sepal.Width  3.4280049 2.871999
#> Petal.Length 1.4620007 4.905992
#> Petal.Width  0.2459998 1.675997
#> 
#> Variances:
#> [,,1]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.15065114  0.13080115   0.02084463  0.01309107
#> Sepal.Width    0.13080115  0.17604529   0.01603245  0.01221458
#> Petal.Length   0.02084463  0.01603245   0.02808260  0.00601568
#> Petal.Width    0.01309107  0.01221458   0.00601568  0.01042365
#> [,,2]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length    0.4000438  0.10865444    0.3994018  0.14368256
#> Sepal.Width     0.1086544  0.10928077    0.1238904  0.07284384
#> Petal.Length    0.3994018  0.12389040    0.6109024  0.25738990
#> Petal.Width     0.1436826  0.07284384    0.2573899  0.16808182

mod2 = densityMclust(faithful, plot = FALSE)
summary(mod2)
#> ------------------------------------------------------- 
#> Density estimation via Gaussian finite mixture modeling 
#> ------------------------------------------------------- 
#> 
#> Mclust EEE (ellipsoidal, equal volume, shape and orientation) model with 3
#> components: 
#> 
#>  log-likelihood   n df       BIC       ICL
#>       -1126.326 272 11 -2314.316 -2357.824
summary(mod2, parameters = TRUE)
#> ------------------------------------------------------- 
#> Density estimation via Gaussian finite mixture modeling 
#> ------------------------------------------------------- 
#> 
#> Mclust EEE (ellipsoidal, equal volume, shape and orientation) model with 3
#> components: 
#> 
#>  log-likelihood   n df       BIC       ICL
#>       -1126.326 272 11 -2314.316 -2357.824
#> 
#> Mixing probabilities:
#>         1         2         3 
#> 0.1656784 0.3563696 0.4779520 
#> 
#> Means:
#>                [,1]      [,2]      [,3]
#> eruptions  3.793066  2.037596  4.463245
#> waiting   77.521051 54.491158 80.833439
#> 
#> Variances:
#> [,,1]
#>            eruptions    waiting
#> eruptions 0.07825448  0.4801979
#> waiting   0.48019785 33.7671464
#> [,,2]
#>            eruptions    waiting
#> eruptions 0.07825448  0.4801979
#> waiting   0.48019785 33.7671464
#> [,,3]
#>            eruptions    waiting
#> eruptions 0.07825448  0.4801979
#> waiting   0.48019785 33.7671464
# }
```
