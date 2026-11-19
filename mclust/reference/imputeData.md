# Missing data imputation via the mix package

Imputes missing data using the mix package.

## Usage

``` r
imputeData(data, categorical = NULL, seed = NULL, verbose = interactive())
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations containing
  missing values. Categorical variables are allowed. If a matrix or data
  frame, rows correspond to observations and columns correspond to
  variables.

- categorical:

  A logical vectors whose *i*th entry is `TRUE` if the *i*th variable or
  column of `data` is to be interpreted as categorical and `FALSE`
  otherwise. The default is to assume that a variable is to be
  interpreted as categorical only if it is a factor.

- seed:

  A seed for the function `rngseed` that is used to initialize the
  random number generator in mix. By default, a seed is chosen uniformly
  in the interval `(.Machine$integer.max/1024, .Machine$integer.max)`.

- verbose:

  A logical, if `TRUE` reports info about iterations of the algorithm.

## Value

A dataset of the same dimensions as `data` with missing values filled
in.

## References

Schafer J. L. (1997). Analysis of Imcomplete Multivariate Data, Chapman
and Hall.

## See also

[`imputePairs`](https://mclust-org.github.io/mclust/reference/imputePairs.md)

## Examples

``` r
# \donttest{
# Note that package 'mix' must be installed
data(stlouis, package = "mix")
 
# impute the continuos variables in the stlouis data
stlimp <- imputeData(stlouis[,-(1:3)])

# plot imputed values
imputePairs(stlouis[,-(1:3)], stlimp)

# }
```
