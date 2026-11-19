# Summarizing semi-supervised classification model based on Gaussian finite mixtures

Summary method for class `"MclustSSC"`.

## Usage

``` r
# S3 method for class 'MclustSSC'
summary(object, parameters = FALSE, ...)
# S3 method for class 'summary.MclustSSC'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- object:

  An object of class `'MclustSSC'` resulting from a call to
  [`MclustSSC`](https://mclust-org.github.io/mclust/reference/MclustSSC.md).

- x:

  An object of class `'summary.MclustSSC'`, usually, a result of a call
  to `summary.MclustSSC`.

- parameters:

  Logical; if `TRUE`, the parameters of mixture components are printed.

- digits:

  The number of significant digits to use when printing.

- ...:

  Further arguments passed to or from other methods.

## Value

The function `summary.MclustSSC` computes and returns a list of summary
statistics of the estimated MclustSSC model for semi-supervised
classification.

## Author

Luca Scrucca

## See also

[`MclustSSC`](https://mclust-org.github.io/mclust/reference/MclustSSC.md),
[`plot.MclustSSC`](https://mclust-org.github.io/mclust/reference/plot.MclustSSC.md).
