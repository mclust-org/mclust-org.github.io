# M-step for parameterized Gaussian mixture models

Maximization step in the EM algorithm for parameterized Gaussian mixture
models.

## Usage

``` r
mstep(data, modelName, z, prior = NULL, warn = NULL, ...)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- modelName:

  A character string indicating the model. The help file for
  [`mclustModelNames`](https://mclust-org.github.io/mclust/reference/mclustModelNames.md)
  describes the available models.

- z:

  A matrix whose `[i,k]`th entry is the conditional probability of the
  ith observation belonging to the *k*th component of the mixture. In
  analyses involving noise, this should not include the conditional
  probabilities for the noise component.

- prior:

  Specification of a conjugate prior on the means and variances. The
  default assumes no prior.

- warn:

  A logical value indicating whether or not certain warnings (usually
  related to singularity) should be issued when the estimation fails.
  The default is given by `mclust.options("warn")`.

- ...:

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A list including the following components:

- modelName:

  A character string identifying the model (same as the input argument).

- parameters:

  `pro`

  : A vector whose *k*th component is the mixing proportion for the
    *k*th component of the mixture model. If the model includes a
    Poisson term for noise, there should be one more mixing proportion
    than the number of Gaussian components.

  `mean`

  : The mean for each component. If there is more than one component,
    this is a matrix whose kth column is the mean of the *k*th component
    of the mixture model.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

- Attributes::

  `"info"` For those models with iterative M-steps (`"VEI"` and
  `"VEV"`), information on the iteration.\
  `"WARNING"` An appropriate warning if problems are encountered in the
  computations.

## Note

This function computes the M-step only for MVN mixtures, so in analyses
involving noise, the conditional probabilities input should exclude
those for the noise component.\

In contrast to `me` for the EM algorithm, computations in `mstep` are
carried out unless failure due to overflow would occur. To impose
stricter tolerances on a single `mstep`, use `me` with the *itmax*
component of the `control` argument set to 1.

## See also

[`mstepE`](https://mclust-org.github.io/mclust/reference/mstepE.md),
...,
[`mstepVVV`](https://mclust-org.github.io/mclust/reference/mstepE.md),
[`emControl`](https://mclust-org.github.io/mclust/reference/emControl.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md).

## Examples

``` r
# \donttest{
mstep(modelName = "VII", data = iris[,-5], z = unmap(iris[,5]))# }
#> $modelName
#> [1] "VII"
#> 
#> $prior
#> NULL
#> 
#> $n
#> [1] 150
#> 
#> $d
#> [1] 4
#> 
#> $G
#> [1] 3
#> 
#> $z
#>        [,1] [,2] [,3]
#>   [1,]    1    0    0
#>   [2,]    1    0    0
#>   [3,]    1    0    0
#>   [4,]    1    0    0
#>   [5,]    1    0    0
#>   [6,]    1    0    0
#>   [7,]    1    0    0
#>   [8,]    1    0    0
#>   [9,]    1    0    0
#>  [10,]    1    0    0
#>  [11,]    1    0    0
#>  [12,]    1    0    0
#>  [13,]    1    0    0
#>  [14,]    1    0    0
#>  [15,]    1    0    0
#>  [16,]    1    0    0
#>  [17,]    1    0    0
#>  [18,]    1    0    0
#>  [19,]    1    0    0
#>  [20,]    1    0    0
#>  [21,]    1    0    0
#>  [22,]    1    0    0
#>  [23,]    1    0    0
#>  [24,]    1    0    0
#>  [25,]    1    0    0
#>  [26,]    1    0    0
#>  [27,]    1    0    0
#>  [28,]    1    0    0
#>  [29,]    1    0    0
#>  [30,]    1    0    0
#>  [31,]    1    0    0
#>  [32,]    1    0    0
#>  [33,]    1    0    0
#>  [34,]    1    0    0
#>  [35,]    1    0    0
#>  [36,]    1    0    0
#>  [37,]    1    0    0
#>  [38,]    1    0    0
#>  [39,]    1    0    0
#>  [40,]    1    0    0
#>  [41,]    1    0    0
#>  [42,]    1    0    0
#>  [43,]    1    0    0
#>  [44,]    1    0    0
#>  [45,]    1    0    0
#>  [46,]    1    0    0
#>  [47,]    1    0    0
#>  [48,]    1    0    0
#>  [49,]    1    0    0
#>  [50,]    1    0    0
#>  [51,]    0    1    0
#>  [52,]    0    1    0
#>  [53,]    0    1    0
#>  [54,]    0    1    0
#>  [55,]    0    1    0
#>  [56,]    0    1    0
#>  [57,]    0    1    0
#>  [58,]    0    1    0
#>  [59,]    0    1    0
#>  [60,]    0    1    0
#>  [61,]    0    1    0
#>  [62,]    0    1    0
#>  [63,]    0    1    0
#>  [64,]    0    1    0
#>  [65,]    0    1    0
#>  [66,]    0    1    0
#>  [67,]    0    1    0
#>  [68,]    0    1    0
#>  [69,]    0    1    0
#>  [70,]    0    1    0
#>  [71,]    0    1    0
#>  [72,]    0    1    0
#>  [73,]    0    1    0
#>  [74,]    0    1    0
#>  [75,]    0    1    0
#>  [76,]    0    1    0
#>  [77,]    0    1    0
#>  [78,]    0    1    0
#>  [79,]    0    1    0
#>  [80,]    0    1    0
#>  [81,]    0    1    0
#>  [82,]    0    1    0
#>  [83,]    0    1    0
#>  [84,]    0    1    0
#>  [85,]    0    1    0
#>  [86,]    0    1    0
#>  [87,]    0    1    0
#>  [88,]    0    1    0
#>  [89,]    0    1    0
#>  [90,]    0    1    0
#>  [91,]    0    1    0
#>  [92,]    0    1    0
#>  [93,]    0    1    0
#>  [94,]    0    1    0
#>  [95,]    0    1    0
#>  [96,]    0    1    0
#>  [97,]    0    1    0
#>  [98,]    0    1    0
#>  [99,]    0    1    0
#> [100,]    0    1    0
#> [101,]    0    0    1
#> [102,]    0    0    1
#> [103,]    0    0    1
#> [104,]    0    0    1
#> [105,]    0    0    1
#> [106,]    0    0    1
#> [107,]    0    0    1
#> [108,]    0    0    1
#> [109,]    0    0    1
#> [110,]    0    0    1
#> [111,]    0    0    1
#> [112,]    0    0    1
#> [113,]    0    0    1
#> [114,]    0    0    1
#> [115,]    0    0    1
#> [116,]    0    0    1
#> [117,]    0    0    1
#> [118,]    0    0    1
#> [119,]    0    0    1
#> [120,]    0    0    1
#> [121,]    0    0    1
#> [122,]    0    0    1
#> [123,]    0    0    1
#> [124,]    0    0    1
#> [125,]    0    0    1
#> [126,]    0    0    1
#> [127,]    0    0    1
#> [128,]    0    0    1
#> [129,]    0    0    1
#> [130,]    0    0    1
#> [131,]    0    0    1
#> [132,]    0    0    1
#> [133,]    0    0    1
#> [134,]    0    0    1
#> [135,]    0    0    1
#> [136,]    0    0    1
#> [137,]    0    0    1
#> [138,]    0    0    1
#> [139,]    0    0    1
#> [140,]    0    0    1
#> [141,]    0    0    1
#> [142,]    0    0    1
#> [143,]    0    0    1
#> [144,]    0    0    1
#> [145,]    0    0    1
#> [146,]    0    0    1
#> [147,]    0    0    1
#> [148,]    0    0    1
#> [149,]    0    0    1
#> [150,]    0    0    1
#> 
#> $parameters
#> $parameters$pro
#> [1] 0.3333333 0.3333333 0.3333333
#> 
#> $parameters$mean
#>               [,1]  [,2]  [,3]
#> Sepal.Length 5.006 5.936 6.588
#> Sepal.Width  3.428 2.770 2.974
#> Petal.Length 1.462 4.260 5.552
#> Petal.Width  0.246 1.326 2.026
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "VII"
#> 
#> $parameters$variance$d
#> [1] 4
#> 
#> $parameters$variance$G
#> [1] 3
#> 
#> $parameters$variance$sigma
#> , , 1
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.075755    0.000000     0.000000    0.000000
#> Sepal.Width      0.000000    0.075755     0.000000    0.000000
#> Petal.Length     0.000000    0.000000     0.075755    0.000000
#> Petal.Width      0.000000    0.000000     0.000000    0.075755
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.153082    0.000000     0.000000    0.000000
#> Sepal.Width      0.000000    0.153082     0.000000    0.000000
#> Petal.Length     0.000000    0.000000     0.153082    0.000000
#> Petal.Width      0.000000    0.000000     0.000000    0.153082
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length      0.21765     0.00000      0.00000     0.00000
#> Sepal.Width       0.00000     0.21765      0.00000     0.00000
#> Petal.Length      0.00000     0.00000      0.21765     0.00000
#> Petal.Width       0.00000     0.00000      0.00000     0.21765
#> 
#> 
#> $parameters$variance$sigmasq
#> [1] 0.075755 0.153082 0.217650
#> 
#> $parameters$variance$scale
#> [1] 0.075755 0.153082 0.217650
#> 
#> 
#> 
#> attr(,"returnCode")
#> [1] 0
```
