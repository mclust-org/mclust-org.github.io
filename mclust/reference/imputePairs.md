# Pairwise Scatter Plots showing Missing Data Imputations

Creates a scatter plot for each pair of variables in given data,
allowing display of imputations for missing values in different colors
and symbols than non missing values.

## Usage

``` r
imputePairs(data, dataImp, 
            symbols = c(1,16), colors = c("black", "red"), labels,
            panel = points, ..., lower.panel = panel, upper.panel = panel, 
            diag.panel = NULL, text.panel = textPanel, label.pos = 0.5 + 
            has.diag/3, cex.labels = NULL, font.labels = 1, row1attop = TRUE, 
            gap = 0.2)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations containing
  missing values. Categorical variables are not allowed. If a matrix or
  data frame, rows correspond to observations and columns correspond to
  variables.

- dataImp:

  The dataset `data` with missing values imputed.

- symbols:

  Either an integer or character vector assigning plotting symbols to
  the nonmissing data and impued values, respectively. The default is a
  closed circle for the nonmissing data and an open circle for the
  imputed values.

- colors:

  Either an integer or character vector assigning colors to the
  nonmissing data and impued values, respectively. The default is black
  for the nonmissing data and red for the imputed values.

- labels:

  As in function `pairs`.

- panel:

  As in function `pairs`.

- ...:

  As in function `pairs`.

- lower.panel:

  As in function `pairs`.

- upper.panel:

  As in function `pairs`.

- diag.panel:

  As in function `pairs`.

- text.panel:

  As in function `pairs`.

- label.pos:

  As in function `pairs`.

- cex.labels:

  As in function `pairs`.

- font.labels:

  As in function `pairs`.

- row1attop:

  As in function `pairs`.

- gap:

  As in function `pairs`.

## Value

A pairs plot displaying the location of missing and nonmissing values.

## References

Schafer J. L. (1997). Analysis of Imcomplete Multivariate Data, Chapman
and Hall.

## See also

[`pairs`](https://rdrr.io/r/graphics/pairs.html),
[`imputeData`](https://mclust-org.github.io/mclust/reference/imputeData.md)

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
