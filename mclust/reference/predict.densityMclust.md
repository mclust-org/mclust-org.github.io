# Density estimate of multivariate observations by Gaussian finite mixture modeling

Compute density estimation for multivariate observations based on
Gaussian finite mixture models estimated by
[`densityMclust`](https://mclust-org.github.io/mclust/reference/densityMclust.md).

## Usage

``` r
# S3 method for class 'densityMclust'
predict(object, newdata, what = c("dens", "cdens", "z"), logarithm = FALSE, ...)
```

## Arguments

- object:

  an object of class `'densityMclust'` resulting from a call to
  [`densityMclust`](https://mclust-org.github.io/mclust/reference/densityMclust.md).

- newdata:

  a vector, a data frame or matrix giving the data. If missing the
  density is computed for the input data obtained from the call to
  [`densityMclust`](https://mclust-org.github.io/mclust/reference/densityMclust.md).

- what:

  a character string specifying what to retrieve: `"dens"` returns a
  vector of values for the mixture density; `"cdens"` returns a matrix
  of component densities for each mixture component (along the columns);
  `"z"` returns a matrix of conditional probabilities of each data point
  to belong to a mixture component.

- logarithm:

  A logical value indicating whether or not the logarithm of the density
  or component densities should be returned.

- ...:

  further arguments passed to or from other methods.

## Value

Returns a vector or a matrix of densities evaluated at `newdata`
depending on the argument `what` (see above).

## Author

Luca Scrucca

## See also

[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md).

## Examples

``` r
# \donttest{
x <- faithful$waiting
dens <- densityMclust(x, plot = FALSE)
x0 <- seq(50, 100, by = 10)
d0 <- predict(dens, x0)
plot(dens, what = "density")
points(x0, d0, pch = 20)

# }
```
