# Classify multivariate observations on a dimension reduced subspace by Gaussian finite mixture modeling

Classify multivariate observations on a dimension reduced subspace
estimated from a Gaussian finite mixture model.

## Usage

``` r
# S3 method for class 'MclustDR'
predict(object, dim = 1:object$numdir, newdata, eval.points, ...)
```

## Arguments

- object:

  an object of class `'MclustDR'` resulting from a call to
  [`MclustDR`](https://mclust-org.github.io/mclust/reference/MclustDR.md)`.`

- dim:

  the dimensions of the reduced subspace used for prediction.

- newdata:

  a data frame or matrix giving the data. If missing the data obtained
  from the call to
  [`MclustDR`](https://mclust-org.github.io/mclust/reference/MclustDR.md)
  are used.

- eval.points:

  a data frame or matrix giving the data projected on the reduced
  subspace. If provided `newdata` is not used.

- ...:

  further arguments passed to or from other methods.

## Value

Returns a list of with the following components:

- dir:

  a matrix containing the data projected onto the `dim` dimensions of
  the reduced subspace.

- density:

  densities from mixture model for each data point.

- z:

  a matrix whose *\[i,k\]*th entry is the probability that observation
  *i* in `newdata` belongs to the *k*th class.

- uncertainty:

  The uncertainty associated with the classification.

- classification:

  A vector of values giving the MAP classification.

## References

Scrucca, L. (2010) Dimension reduction for model-based clustering.
*Statistics and Computing*, 20(4), pp. 471-484.

## Author

Luca Scrucca

## See also

[`MclustDR`](https://mclust-org.github.io/mclust/reference/MclustDR.md).

## Examples

``` r
mod = Mclust(iris[,1:4])
dr = MclustDR(mod)
pred = predict(dr)
str(pred)
#> List of 5
#>  $ dir           : num [1:150, 1] 1.89 1.49 1.66 1.44 1.94 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr "Dir1"
#>  $ density       : num [1:150] 0.363 0.28 0.366 0.245 0.341 ...
#>  $ z             : num [1:150, 1:2] 1 1 1 1 1 ...
#>  $ uncertainty   : num [1:150] 1.08e-09 3.75e-07 2.88e-08 8.46e-07 5.17e-10 ...
#>  $ classification: Factor w/ 2 levels "1","2": 1 1 1 1 1 1 1 1 1 1 ...

data(banknote)
mod = MclustDA(banknote[,2:7], banknote$Status)
dr = MclustDR(mod)
pred = predict(dr)
str(pred)
#> List of 5
#>  $ dir           : num [1:200, 1:2] -0.946 -1.751 -1.583 -1.742 -1.528 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "Dir1" "Dir2"
#>  $ density       : num [1:200] 0.137 0.277 0.169 0.31 0.462 ...
#>  $ z             : num [1:200, 1:2] 1.47e-08 4.43e-16 4.83e-18 2.49e-17 1.67e-15 ...
#>  $ uncertainty   : num [1:200] 1.47e-08 4.44e-16 0.00 0.00 1.78e-15 ...
#>  $ classification: Factor w/ 2 levels "counterfeit",..: 2 2 2 2 2 2 2 2 2 2 ...
```
