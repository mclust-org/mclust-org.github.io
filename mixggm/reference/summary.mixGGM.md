# Summarizing Mixture of Gaussian Graphical Models

Summary method for class `'mixGGM'`.

## Usage

``` r
# S3 method for class 'mixGGM'
summary(object, graphs = TRUE, clusters = FALSE, parameters = FALSE, ...)
# S3 method for class 'summary.mixGGM'
print(x, digits = getOption("digits"), ...)
```

## Arguments

- object:

  An object of class `'mixGGM'` resulting of a call to
  [`mixGGM`](https://mclust-org.github.io/mixggm/reference/mixGGM.md).

- x:

  An object of class `'summary.mixGGM'`, usually, a result of a call to
  `summary.mixGGM`.

- graphs:

  Logical; if `TRUE`, the adjancecy matrices of mixture graphs is
  printed.

- clusters:

  Logical; if `TRUE`, the table of MAP classification/clustering of
  observations is printed.

- parameters:

  Logical; if `TRUE`, the parameters of mixture components are printed.

- digits:

  The number of significant digits to use when printing.

- ...:

  Further arguments passed to or from other methods.

## See also

[`mixGGM`](https://mclust-org.github.io/mixggm/reference/mixGGM.md).
