# Summarizing dimension reduction method for model-based clustering and classification

Summary method for class `"MclustDR"`.

## Usage

``` r
# S3 method for class 'MclustDR'
summary(object, numdir, std = FALSE, ...)
# S3 method for class 'summary.MclustDR'
print(x, digits = max(5, getOption("digits") - 3), ...)
```

## Arguments

- object:

  An object of class `'MclustDR'` resulting from a call to
  [`MclustDR`](https://mclust-org.github.io/mclust/reference/MclustDR.md).

- x:

  An object of class `'summary.MclustDR'`, usually, a result of a call
  to `summary.MclustDR`.

- numdir:

  An integer providing the number of basis directions to be printed.

- std:

  if `TRUE` the coefficients basis are scaled such that all predictors
  have unit standard deviation.

- digits:

  The number of significant digits to use when printing.

- ...:

  Further arguments passed to or from other methods.

## Author

Luca Scrucca

## See also

[`MclustDR`](https://mclust-org.github.io/mclust/reference/MclustDR.md),
[`plot.MclustDR`](https://mclust-org.github.io/mclust/reference/plot.MclustDR.md)
