# EM algorithm starting with E-step for a parameterized Gaussian mixture model

Implements the EM algorithm for a parameterized Gaussian mixture model,
starting with the expectation step.

## Usage

``` r
emE(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emV(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emX(data, prior = NULL, warn = NULL, ...)
emEII(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVII(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emEEI(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVEI(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emEVI(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVVI(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emEEE(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVEE(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emEVE(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVVE(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emEEV(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVEV(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emEVV(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emVVV(data, parameters, prior = NULL, control = emControl(), warn = NULL, ...)
emXII(data, prior = NULL, warn = NULL, ...)
emXXI(data, prior = NULL, warn = NULL, ...)
emXXX(data, prior = NULL, warn = NULL, ...)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- parameters:

  The parameters of the model:

  `pro`

  : Mixing proportions for the components of the mixture. There should
    one more mixing proportion than the number of Gaussian components if
    the mixture model includes a Poisson noise term.

  `mean`

  : The mean for each component. If there is more than one component,
    this is a matrix whose kth column is the mean of the *k*th component
    of the mixture model.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

  `Vinv`

  : An estimate of the reciprocal hypervolume of the data region. The
    default is determined by applying function `hypvol` to the data.
    Used only when `pro` includes an additional mixing proportion for a
    noise component.

- prior:

  The default assumes no prior, but this argument allows specification
  of a conjugate prior on the means and variances through the function
  `priorControl`.

- control:

  A list of control parameters for EM. The defaults are set by the call
  [`emControl()`](https://mclust-org.github.io/mclust/reference/emControl.md).

- warn:

  A logical value indicating whether or not a warning should be issued
  whenever a singularity is encountered. The default is given in
  `mclust.options("warn")`.

- ...:

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A list including the following components:

- modelName:

  A character string identifying the model (same as the input argument).

- z:

  A matrix whose `[i,k]`th entry is the conditional probability of the
  *i*th observation belonging to the *k*th component of the mixture.

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

  `Vinv`

  : The estimate of the reciprocal hypervolume of the data region used
    in the computation when the input indicates the addition of a noise
    component to the model.

- loglik:

  The log likelihood for the data in the mixture model.

- Attributes::

  `"info"` Information on the iteration.\
  `"WARNING"` An appropriate warning if problems are encountered in the
  computations.

## See also

[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md).

## Examples

``` r
# \donttest{
msEst <- mstepEEE(data = iris[,-5], z = unmap(iris[,5]))
names(msEst)
#> [1] "modelName"  "prior"      "n"          "d"          "G"         
#> [6] "z"          "parameters"

emEEE(data = iris[,-5], parameters = msEst$parameters)# }
#> $modelName
#> [1] "EEE"
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
#>                [,1]         [,2]         [,3]
#>   [1,] 1.000000e+00 5.669389e-22 2.970668e-42
#>   [2,] 1.000000e+00 7.125511e-18 2.890550e-37
#>   [3,] 1.000000e+00 2.209921e-19 5.645717e-39
#>   [4,] 1.000000e+00 1.944369e-16 4.493925e-35
#>   [5,] 1.000000e+00 2.969620e-22 1.819554e-42
#>   [6,] 1.000000e+00 7.746074e-21 1.005427e-39
#>   [7,] 1.000000e+00 2.519429e-18 5.970904e-37
#>   [8,] 1.000000e+00 5.638578e-20 1.231056e-39
#>   [9,] 1.000000e+00 2.778048e-15 1.088991e-33
#>  [10,] 1.000000e+00 1.193692e-18 1.753502e-38
#>  [11,] 1.000000e+00 1.635161e-23 3.372329e-44
#>  [12,] 1.000000e+00 2.937431e-18 3.124729e-37
#>  [13,] 1.000000e+00 1.517011e-18 1.969537e-38
#>  [14,] 1.000000e+00 1.843982e-19 1.764184e-39
#>  [15,] 1.000000e+00 6.849292e-30 1.303566e-52
#>  [16,] 1.000000e+00 3.245264e-27 7.704284e-48
#>  [17,] 1.000000e+00 1.264758e-24 7.386154e-45
#>  [18,] 1.000000e+00 6.462169e-21 1.584880e-40
#>  [19,] 1.000000e+00 2.203157e-22 2.328338e-42
#>  [20,] 1.000000e+00 5.749059e-22 1.456245e-41
#>  [21,] 1.000000e+00 1.272471e-19 2.601003e-39
#>  [22,] 1.000000e+00 3.035910e-20 4.610963e-39
#>  [23,] 1.000000e+00 1.681457e-24 2.334188e-45
#>  [24,] 1.000000e+00 1.247580e-14 1.126055e-31
#>  [25,] 1.000000e+00 2.033636e-15 2.214427e-33
#>  [26,] 1.000000e+00 2.297938e-16 2.933735e-35
#>  [27,] 1.000000e+00 6.480706e-17 6.730469e-35
#>  [28,] 1.000000e+00 2.066760e-21 1.569681e-41
#>  [29,] 1.000000e+00 1.082360e-21 4.850018e-42
#>  [30,] 1.000000e+00 1.529968e-16 4.000996e-35
#>  [31,] 1.000000e+00 2.920907e-16 6.532168e-35
#>  [32,] 1.000000e+00 2.112491e-19 2.006593e-38
#>  [33,] 1.000000e+00 1.833782e-26 6.732546e-48
#>  [34,] 1.000000e+00 3.568839e-28 6.559131e-50
#>  [35,] 1.000000e+00 1.360613e-17 9.355103e-37
#>  [36,] 1.000000e+00 1.748091e-21 6.118695e-42
#>  [37,] 1.000000e+00 1.848029e-24 8.856622e-46
#>  [38,] 1.000000e+00 6.322278e-23 1.239791e-43
#>  [39,] 1.000000e+00 6.778306e-17 9.552703e-36
#>  [40,] 1.000000e+00 2.323565e-20 3.386508e-40
#>  [41,] 1.000000e+00 1.772656e-21 2.999433e-41
#>  [42,] 1.000000e+00 1.458475e-11 3.636199e-29
#>  [43,] 1.000000e+00 3.158065e-18 2.712042e-37
#>  [44,] 1.000000e+00 1.817425e-15 3.227861e-32
#>  [45,] 1.000000e+00 4.013405e-17 1.057569e-34
#>  [46,] 1.000000e+00 1.970937e-16 5.605941e-35
#>  [47,] 1.000000e+00 4.461935e-22 5.242944e-42
#>  [48,] 1.000000e+00 4.744169e-18 3.942102e-37
#>  [49,] 1.000000e+00 3.968032e-23 1.225902e-43
#>  [50,] 1.000000e+00 2.952916e-20 3.803732e-40
#>  [51,] 4.352061e-18 9.999271e-01 7.285356e-05
#>  [52,] 1.462104e-19 9.991833e-01 8.167043e-04
#>  [53,] 4.326581e-22 9.969243e-01 3.075740e-03
#>  [54,] 3.588876e-22 9.997272e-01 2.727797e-04
#>  [55,] 8.685517e-23 9.968195e-01 3.180529e-03
#>  [56,] 8.313347e-23 9.983021e-01 1.697924e-03
#>  [57,] 3.067079e-22 9.793448e-01 2.065525e-02
#>  [58,] 5.113810e-14 9.999999e-01 1.031644e-07
#>  [59,] 1.272609e-19 9.999241e-01 7.588810e-05
#>  [60,] 8.975128e-21 9.992664e-01 7.336194e-04
#>  [61,] 3.442085e-18 9.999991e-01 8.743942e-07
#>  [62,] 5.603197e-20 9.990125e-01 9.874887e-04
#>  [63,] 9.656653e-18 9.999995e-01 4.518974e-07
#>  [64,] 1.491160e-23 9.942156e-01 5.784381e-03
#>  [65,] 5.275521e-14 9.999981e-01 1.854263e-06
#>  [66,] 4.551187e-17 9.999694e-01 3.064952e-05
#>  [67,] 5.482635e-24 9.671064e-01 3.289363e-02
#>  [68,] 3.956205e-16 9.999994e-01 6.381923e-07
#>  [69,] 5.338555e-27 9.786346e-01 2.136543e-02
#>  [70,] 2.149035e-17 9.999977e-01 2.333318e-06
#>  [71,] 2.269907e-28 1.338383e-01 8.661617e-01
#>  [72,] 1.564633e-16 9.999930e-01 6.999013e-06
#>  [73,] 1.862337e-28 8.664919e-01 1.335081e-01
#>  [74,] 4.204653e-22 9.996599e-01 3.400895e-04
#>  [75,] 1.496225e-17 9.999834e-01 1.663712e-05
#>  [76,] 4.048066e-18 9.999412e-01 5.881490e-05
#>  [77,] 1.811813e-22 9.990449e-01 9.551325e-04
#>  [78,] 9.733911e-27 7.032941e-01 2.967059e-01
#>  [79,] 4.207095e-23 9.914214e-01 8.578608e-03
#>  [80,] 1.686652e-11 1.000000e+00 1.168802e-08
#>  [81,] 1.691014e-17 9.999979e-01 2.062222e-06
#>  [82,] 1.705135e-15 9.999998e-01 2.029185e-07
#>  [83,] 2.383037e-16 9.999970e-01 2.965664e-06
#>  [84,] 4.543859e-32 1.255609e-01 8.744391e-01
#>  [85,] 8.944225e-25 9.290885e-01 7.091151e-02
#>  [86,] 7.853508e-21 9.883971e-01 1.160287e-02
#>  [87,] 5.759078e-21 9.985336e-01 1.466369e-03
#>  [88,] 7.047203e-23 9.997609e-01 2.391093e-04
#>  [89,] 4.510692e-18 9.999302e-01 6.984723e-05
#>  [90,] 7.703797e-21 9.998338e-01 1.662373e-04
#>  [91,] 6.639378e-23 9.993841e-01 6.159223e-04
#>  [92,] 6.134157e-22 9.979127e-01 2.087322e-03
#>  [93,] 5.814465e-18 9.999918e-01 8.248956e-06
#>  [94,] 2.678609e-14 9.999999e-01 8.822342e-08
#>  [95,] 5.126445e-21 9.996812e-01 3.187501e-04
#>  [96,] 1.410433e-17 9.999784e-01 2.163084e-05
#>  [97,] 2.670627e-19 9.998703e-01 1.296855e-04
#>  [98,] 2.540730e-18 9.999627e-01 3.733285e-05
#>  [99,] 8.473985e-11 1.000000e+00 1.640967e-08
#> [100,] 5.099820e-19 9.999235e-01 7.651826e-05
#> [101,] 8.752565e-53 1.844658e-09 1.000000e+00
#> [102,] 2.587892e-38 6.236363e-04 9.993764e-01
#> [103,] 1.986727e-42 2.318538e-05 9.999768e-01
#> [104,] 1.180532e-38 7.486688e-04 9.992513e-01
#> [105,] 3.099804e-46 9.518578e-07 9.999990e-01
#> [106,] 1.307361e-48 7.687368e-07 9.999992e-01
#> [107,] 9.267670e-34 2.246431e-02 9.775357e-01
#> [108,] 4.933987e-42 1.873452e-04 9.998127e-01
#> [109,] 4.505583e-42 2.970581e-04 9.997029e-01
#> [110,] 1.055915e-46 6.784148e-08 9.999999e-01
#> [111,] 2.974828e-32 7.725660e-03 9.922743e-01
#> [112,] 1.617173e-37 1.493449e-03 9.985066e-01
#> [113,] 5.629068e-39 1.532789e-04 9.998467e-01
#> [114,] 7.280330e-41 1.177402e-04 9.998823e-01
#> [115,] 3.555665e-46 3.558495e-07 9.999996e-01
#> [116,] 1.471997e-40 1.075125e-05 9.999892e-01
#> [117,] 1.771439e-35 4.654684e-03 9.953453e-01
#> [118,] 7.138648e-45 8.218517e-07 9.999992e-01
#> [119,] 1.899015e-58 1.906714e-09 1.000000e+00
#> [120,] 5.082160e-33 2.972470e-01 7.027530e-01
#> [121,] 6.864466e-43 3.648876e-06 9.999964e-01
#> [122,] 8.039862e-38 3.586891e-04 9.996413e-01
#> [123,] 3.747574e-49 1.512670e-06 9.999985e-01
#> [124,] 2.931140e-31 9.409549e-02 9.059045e-01
#> [125,] 8.774784e-40 4.563353e-05 9.999544e-01
#> [126,] 2.005765e-36 2.684624e-03 9.973154e-01
#> [127,] 8.506186e-30 1.616826e-01 8.383174e-01
#> [128,] 4.664639e-30 8.868082e-02 9.113192e-01
#> [129,] 4.765225e-44 8.543916e-06 9.999915e-01
#> [130,] 5.127429e-32 1.448711e-01 8.551289e-01
#> [131,] 2.089865e-41 2.206495e-04 9.997794e-01
#> [132,] 1.902049e-36 4.134090e-04 9.995866e-01
#> [133,] 8.931905e-46 1.825407e-06 9.999982e-01
#> [134,] 2.029625e-28 7.432053e-01 2.567947e-01
#> [135,] 3.219380e-35 7.082646e-02 9.291735e-01
#> [136,] 4.365673e-45 2.536647e-06 9.999975e-01
#> [137,] 3.772459e-45 2.458322e-07 9.999998e-01
#> [138,] 2.894065e-35 3.983245e-03 9.960168e-01
#> [139,] 2.370290e-29 1.236116e-01 8.763884e-01
#> [140,] 2.331576e-36 6.383548e-04 9.993616e-01
#> [141,] 3.151247e-45 5.888272e-07 9.999994e-01
#> [142,] 5.807134e-36 2.983685e-04 9.997016e-01
#> [143,] 2.587892e-38 6.236363e-04 9.993764e-01
#> [144,] 5.118189e-46 5.166787e-07 9.999995e-01
#> [145,] 1.083147e-46 9.508309e-08 9.999999e-01
#> [146,] 3.855859e-39 4.781458e-05 9.999522e-01
#> [147,] 8.908459e-36 6.189358e-03 9.938106e-01
#> [148,] 4.421495e-35 2.180277e-03 9.978197e-01
#> [149,] 2.042704e-41 3.621223e-06 9.999964e-01
#> [150,] 1.040288e-33 9.114476e-03 9.908855e-01
#> 
#> $parameters
#> $parameters$pro
#> [1] 0.3333333 0.3295291 0.3371375
#> 
#> $parameters$mean
#>               [,1]     [,2]     [,3]
#> Sepal.Length 5.006 5.942204 6.574579
#> Sepal.Width  3.428 2.760815 2.980676
#> Petal.Length 1.462 4.258445 5.538942
#> Petal.Width  0.246 1.319151 2.024796
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "EEE"
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
#> Sepal.Length   0.26391759  0.08986519   0.16961981  0.03933757
#> Sepal.Width    0.08986519  0.11196054   0.05115079  0.03000199
#> Petal.Length   0.16961981  0.05115079   0.18645154  0.04196883
#> Petal.Width    0.03933757  0.03000199   0.04196883  0.03973224
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.26391759  0.08986519   0.16961981  0.03933757
#> Sepal.Width    0.08986519  0.11196054   0.05115079  0.03000199
#> Petal.Length   0.16961981  0.05115079   0.18645154  0.04196883
#> Petal.Width    0.03933757  0.03000199   0.04196883  0.03973224
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.26391759  0.08986519   0.16961981  0.03933757
#> Sepal.Width    0.08986519  0.11196054   0.05115079  0.03000199
#> Petal.Length   0.16961981  0.05115079   0.18645154  0.04196883
#> Petal.Width    0.03933757  0.03000199   0.04196883  0.03973224
#> 
#> 
#> $parameters$variance$Sigma
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.26391759  0.08986519   0.16961981  0.03933757
#> Sepal.Width    0.08986519  0.11196054   0.05115079  0.03000199
#> Petal.Length   0.16961981  0.05115079   0.18645154  0.04196883
#> Petal.Width    0.03933757  0.03000199   0.04196883  0.03973224
#> 
#> $parameters$variance$cholSigma
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   -0.5137291  -0.1749272  -0.33017364 -0.07657259
#> Sepal.Width     0.0000000   0.2852385  -0.02315803  0.05822273
#> Petal.Length    0.0000000   0.0000000  -0.27730960 -0.06503526
#> Petal.Width     0.0000000   0.0000000   0.00000000  0.16201670
#> 
#> 
#> 
#> attr(,"returnCode")
#> [1] 0
```
