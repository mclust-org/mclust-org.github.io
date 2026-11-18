# Cluster prediction by Mixture of Gaussian Graphical Models

Cluster prediction for multivariate observations based on Mixture of
Gaussian Graphical Models estimated by [`mixGGM`](mixGGM.md).

## Usage

``` r
# S3 method for class 'mixGGM'
predict(object, newdata, ...)
```

## Arguments

- object:

  An object of class `'mixGGM'` resulting from a call to
  [`mixGGM`](mixGGM.md).

- newdata:

  A data frame or matrix giving the data. If missing the clustering data
  obtained from the call to [`mixGGM`](mixGGM.md) are classified.

- ...:

  Further arguments passed to or from other methods.

## Value

Returns a list of with the following components:

- classification:

  a factor of predicted cluster labels for `newdata`.

- z:

  a matrix whose *\[i,k\]*th entry is the probability that observation
  *i* in `newdata` belongs to the *k*th cluster.

## See also

[`mixGGM`](mixGGM.md).
