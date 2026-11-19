# Deprecated Functions in mclust package

These functions are provided for compatibility with older versions of
the mclust package only, and may be removed eventually.

## Usage

``` r
cv.MclustDA(...)
cv1EMtrain(data, labels, modelNames=NULL)
bicEMtrain(data, labels, modelNames=NULL)
```

## Arguments

- ...:

  pass arguments down.

- data:

  A numeric vector or matrix of observations.

- labels:

  Labels for each element or row in the dataset.

- modelNames:

  Vector of model names that should be tested. The default is to select
  all available model names.

## See also

[`deprecated`](https://rdrr.io/r/base/Deprecated.html)
