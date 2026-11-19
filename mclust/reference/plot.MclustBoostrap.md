# Plot of bootstrap distributions for mixture model parameters

Plots the bootstrap distribution of parameters as returned by the
[`MclustBootstrap`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md)
function.

## Usage

``` r
# S3 method for class 'MclustBootstrap'
plot(x, what = c("pro", "mean", "var"), 
     show.parest = TRUE, show.confint = TRUE,
     hist.col = "grey", hist.border = "lightgrey", breaks = NA, 
     col = "forestgreen", lwd = 2, lty = 3, 
     xlab = NULL, xlim = NULL, ylim = NULL, ...)
```

## Arguments

- x:

  Object returned by `MclustBootstrap`.

- what:

  Character string specifying if mixing proportions (`"pro"`), component
  means (`"mean"`) or component variances (`"var"`) should be drawn.

- show.parest:

  A logical specifying if the parameter estimate should be drawn as
  vertical line.

- show.confint:

  A logical specifying if the resampling-based confidence interval
  should be drawn at the bottom of the graph. Confidence level can be
  provided as further argument `conf.level`; see
  [`summary.MclustBootstrap`](https://mclust-org.github.io/mclust/reference/summary.MclustBootstrap.md).

- hist.col:

  The color to be used to fill the bars of the histograms.

- hist.border:

  The color of the border around the bars of the histograms.

- breaks:

  The number of breaks used in histogram to visualize the bootstrap
  distribution. When `NA` (default), it is computed as the maximum
  between the "Sturges" and "Freedman-Diaconis" methods; see
  [`nclass.Sturges`](https://rdrr.io/r/grDevices/nclass.html) and
  [`nclass.FD`](https://rdrr.io/r/grDevices/nclass.html).

- col, lwd, lty:

  The color, line width and line type to be used to represent the
  estimated parameters and confidence intervals.

- xlab:

  Optional label for the horizontal axis.

- xlim, ylim:

  A two-values vector of axis range for, respectively, horizontal and
  vertical axis.

- ...:

  Other graphics parameters.

## Value

A plot for each variable/component of the selected parameters.

## See also

[`MclustBootstrap`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md)

## Examples

``` r
# \donttest{
data(diabetes)
X <- diabetes[,-1]
modClust <- Mclust(X, G = 3, modelNames = "VVV")
bootClust <- MclustBootstrap(modClust, nboot = 99)
par(mfrow = c(1,3), mar = c(4,2,2,0.5))
plot(bootClust, what = "pro")

par(mfrow = c(3,3), mar = c(4,2,2,0.5))
plot(bootClust, what = "mean")

# }
```
