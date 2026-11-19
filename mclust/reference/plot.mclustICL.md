# ICL Plot for Model-Based Clustering

Plots the ICL values returned by the
[`mclustICL`](https://mclust-org.github.io/mclust/reference/mclustICL.md)
function.

## Usage

``` r
# S3 method for class 'mclustICL'
plot(x, ylab = "ICL", ...)
```

## Arguments

- x:

  Output from
  [`mclustICL`](https://mclust-org.github.io/mclust/reference/mclustICL.md).

- ylab:

  Label for the vertical axis of the plot.

- ...:

  Further arguments passed to the
  [`plot.mclustBIC`](https://mclust-org.github.io/mclust/reference/plot.mclustBIC.md)
  function.

## Value

A plot of the ICL values.

## See also

[`mclustICL`](https://mclust-org.github.io/mclust/reference/mclustICL.md)

## Examples

``` r
# \donttest{
data(faithful)
faithful.ICL = mclustICL(faithful)
plot(faithful.ICL)

# }
```
