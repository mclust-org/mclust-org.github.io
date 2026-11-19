# Plotting method for Mclust model-based clustering

Plots for model-based clustering results, such as BIC, classification,
uncertainty and density.

## Usage

``` r
# S3 method for class 'Mclust'
plot(x, what = c("BIC", "classification", "uncertainty", "density"), 
     dimens = NULL, xlab = NULL, ylab = NULL,
     addEllipses = TRUE, main = FALSE, ...)
```

## Arguments

- x:

  Output from `Mclust`.

- what:

  A string specifying the type of graph requested. Available choices
  are:

  `"BIC"`

  : plot of BIC values used for choosing the number of clusters.

  `"classification"` =

  : a plot showing the clustering. For data in more than two dimensions
    a pairs plot is produced, followed by a coordinate projection plot
    using specified `dimens`. Ellipses corresponding to covariances of
    mixture components are also drawn if `addEllipses = TRUE`.

  `"uncertainty"`

  : a plot of classification uncertainty. For data in more than two
    dimensions a coordinate projection plot is drawn using specified
    `dimens`.

  `"density"`

  : a plot of estimated density. For data in more than two dimensions a
    matrix of contours for coordinate projection plot is drawn using
    specified `dimens`.

  If not specified, in interactive sessions a menu of choices is
  proposed.

- dimens:

  A vector of integers specifying the dimensions of the coordinate
  projections in case of `"classification"`, `"uncertainty"`, or
  `"density"` plots.

- xlab, ylab:

  Optional labels for the x-axis and the y-axis.

- addEllipses:

  A logical indicating whether or not to add ellipses with axes
  corresponding to the within-cluster covariances in case of
  `"classification"` or `"uncertainty"` plots.

- main:

  A logical or `NULL` indicating whether or not to add a title to the
  plot identifying the type of plot drawn.

- ...:

  Other graphics parameters.

## Details

For more flexibility in plotting, use `mclust1Dplot`, `mclust2Dplot`,
`surfacePlot`, `coordProj`, or `randProj`.

## See also

[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md),
[`plot.mclustBIC`](https://mclust-org.github.io/mclust/reference/plot.mclustBIC.md),
[`plot.mclustICL`](https://mclust-org.github.io/mclust/reference/plot.mclustICL.md),
[`mclust1Dplot`](https://mclust-org.github.io/mclust/reference/mclust1Dplot.md),
[`mclust2Dplot`](https://mclust-org.github.io/mclust/reference/mclust2Dplot.md),
[`surfacePlot`](https://mclust-org.github.io/mclust/reference/surfacePlot.md),
[`coordProj`](https://mclust-org.github.io/mclust/reference/coordProj.md),
[`randProj`](https://mclust-org.github.io/mclust/reference/randProj.md).

## Examples

``` r
# \donttest{
precipMclust <- Mclust(precip)
plot(precipMclust)





faithfulMclust <- Mclust(faithful)
plot(faithfulMclust)





irisMclust <- Mclust(iris[,-5])
plot(irisMclust)




# }
```
