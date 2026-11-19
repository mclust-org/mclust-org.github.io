# Summarizing discriminant analysis based on Gaussian finite mixture modeling

Summary method for class `"MclustDA"`.

## Usage

``` r
# S3 method for class 'MclustDA'
summary(object, parameters = FALSE, newdata, newclass, ...)
# S3 method for class 'summary.MclustDA'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- object:

  An object of class `'MclustDA'` resulting from a call to
  [`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

- x:

  An object of class `'summary.MclustDA'`, usually, a result of a call
  to `summary.MclustDA`.

- parameters:

  Logical; if `TRUE`, the parameters of mixture components are printed.

- newdata:

  A data frame or matrix giving the test data.

- newclass:

  A vector giving the class labels for the observations in the test
  data.

- digits:

  The number of significant digits to use when printing.

- ...:

  Further arguments passed to or from other methods.

## Value

The function `summary.MclustDA` computes and returns a list of summary
statistics of the estimated MclustDA or EDDA model for classification.

## Author

Luca Scrucca

## See also

[`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md),
[`plot.MclustDA`](https://mclust-org.github.io/mclust/reference/plot.MclustDA.md).

## Examples

``` r
mod = MclustDA(data = iris[,1:4], class = iris$Species)
summary(mod)
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
summary(mod, parameters = TRUE)
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
#> Parameter estimates
#> ------------------------------------------------ 
#> Class prior probabilities:
#>     setosa versicolor  virginica 
#>  0.3333333  0.3333333  0.3333333 
#> 
#> Class = setosa
#> -------------- 
#> 
#> Mixing probabilities: 0.8081091 0.1918909 
#> 
#> Means:
#>                   [,1]      [,2]
#> Sepal.Length 4.9484622 5.2483085
#> Sepal.Width  3.3627745 3.7026840
#> Petal.Length 1.4322810 1.5871556
#> Petal.Width  0.2036194 0.4244774
#> 
#> Variances:
#> [,,1]
#>               Sepal.Length  Sepal.Width Petal.Length   Petal.Width
#> Sepal.Length  0.1078221133  0.081427298 0.0088268249 -0.0001451869
#> Sepal.Width   0.0814272976  0.122899586 0.0033006551 -0.0025292826
#> Petal.Length  0.0088268249  0.003300655 0.0258364946  0.0006438246
#> Petal.Width  -0.0001451869 -0.002529283 0.0006438246  0.0033200166
#> [,,2]
#>               Sepal.Length  Sepal.Width Petal.Length   Petal.Width
#> Sepal.Length  0.1078221133  0.081427298 0.0088268249 -0.0001451869
#> Sepal.Width   0.0814272976  0.122899586 0.0033006551 -0.0025292826
#> Petal.Length  0.0088268249  0.003300655 0.0258364946  0.0006438246
#> Petal.Width  -0.0001451869 -0.002529283 0.0006438246  0.0033200166
#> 
#> Class = versicolor
#> ------------------ 
#> 
#> Mixing probabilities: 1 
#> 
#> Means:
#>               [,1]
#> Sepal.Length 5.936
#> Sepal.Width  2.770
#> Petal.Length 4.260
#> Petal.Width  1.326
#> 
#> Variances:
#> [,,1]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.261104     0.08348      0.17924    0.054664
#> Sepal.Width      0.083480     0.09650      0.08100    0.040380
#> Petal.Length     0.179240     0.08100      0.21640    0.071640
#> Petal.Width      0.054664     0.04038      0.07164    0.038324
#> 
#> Class = virginica
#> ----------------- 
#> 
#> Mixing probabilities: 1 
#> 
#> Means:
#>               [,1]
#> Sepal.Length 6.588
#> Sepal.Width  2.974
#> Petal.Length 5.552
#> Petal.Width  2.026
#> 
#> Variances:
#> [,,1]
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.396256    0.091888     0.297224    0.048112
#> Sepal.Width      0.091888    0.101924     0.069952    0.046676
#> Petal.Length     0.297224    0.069952     0.298496    0.047848
#> Petal.Width      0.048112    0.046676     0.047848    0.073924
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
```
