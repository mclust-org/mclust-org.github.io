# E-step in the EM algorithm for a parameterized Gaussian mixture model.

Implements the expectation step in the EM algorithm for a parameterized
Gaussian mixture model.

## Usage

``` r
estepE(data, parameters, warn = NULL, ...)
estepV(data, parameters, warn = NULL, ...)
estepEII(data, parameters, warn = NULL, ...)
estepVII(data, parameters, warn = NULL, ...)
estepEEI(data, parameters, warn = NULL, ...)
estepVEI(data, parameters, warn = NULL, ...)
estepEVI(data, parameters, warn = NULL, ...)
estepVVI(data, parameters, warn = NULL, ...)
estepEEE(data, parameters, warn = NULL, ...)
estepEEV(data, parameters, warn = NULL, ...)
estepVEV(data, parameters, warn = NULL, ...)
estepVVV(data, parameters, warn = NULL, ...)
estepEVE(data, parameters, warn = NULL, ...)
estepEVV(data, parameters, warn = NULL, ...)
estepVEE(data, parameters, warn = NULL, ...)
estepVVE(data, parameters, warn = NULL, ...)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- parameters:

  The parameters of the model:

  `pro`

  : Mixing proportions for the components of the mixture. If the model
    includes a Poisson term for noise, there should be one more mixing
    proportion than the number of Gaussian components.

  mu

  : The mean for each component. If there is more than one component,
    this is a matrix whose columns are the means of the components.

  `variance`

  : A list of variance parameters for the model. The components of this
    list depend on the model specification. See the help file for
    [`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)
    for details.

  `Vinv`

  : An estimate of the reciprocal hypervolume of the data region. If not
    supplied or set to a negative value, the default is determined by
    applying function `hypvol` to the data. Used only when `pro`
    includes an additional mixing proportion for a noise component.

- warn:

  A logical value indicating whether or certain warnings should be
  issued. The default is given by `mclust.options("warn")`.

- ...:

  Catches unused arguments in indirect or list calls via `do.call`.

## Value

A list including the following components:

- modelName:

  Character string identifying the model.

- z:

  A matrix whose `[i,k]`th entry is the conditional probability of the
  *i*th observation belonging to the *k*th component of the mixture.

- parameters:

  The input parameters.

- loglik:

  The logliklihood for the data in the mixture model.

- Attribute:

  `"WARNING"`: An appropriate warning if problems are encountered in the
  computations.

## See also

[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`em`](https://mclust-org.github.io/mclust/reference/em.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`do.call`](https://rdrr.io/r/base/do.call.html),
[`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md).

## Examples

``` r
# \donttest{
msEst <- mstepEII(data = iris[,-5], z = unmap(iris[,5]))
names(msEst)
#> [1] "modelName"  "prior"      "n"          "d"          "G"         
#> [6] "z"          "parameters"

estepEII(data = iris[,-5], parameters = msEst$parameters)# }
#> $modelName
#> [1] "EII"
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
#>   [1,] 1.000000e+00 2.803384e-16 2.385590e-34
#>   [2,] 1.000000e+00 7.327699e-16 1.308337e-34
#>   [3,] 1.000000e+00 1.323453e-17 5.432319e-37
#>   [4,] 1.000000e+00 4.734701e-16 6.206020e-35
#>   [5,] 1.000000e+00 9.644772e-17 6.074008e-35
#>   [6,] 1.000000e+00 3.745645e-13 7.110633e-29
#>   [7,] 1.000000e+00 3.962220e-17 5.263751e-36
#>   [8,] 1.000000e+00 1.530343e-15 1.745542e-33
#>   [9,] 1.000000e+00 5.012695e-17 8.729744e-37
#>  [10,] 1.000000e+00 1.493770e-15 4.553180e-34
#>  [11,] 1.000000e+00 4.946353e-15 4.909613e-32
#>  [12,] 1.000000e+00 2.874117e-15 3.251952e-33
#>  [13,] 1.000000e+00 1.898603e-16 1.366666e-35
#>  [14,] 1.000000e+00 2.965476e-20 1.765968e-41
#>  [15,] 1.000000e+00 5.679989e-17 3.628138e-34
#>  [16,] 1.000000e+00 6.232567e-15 1.539704e-30
#>  [17,] 1.000000e+00 2.030556e-16 1.196559e-33
#>  [18,] 1.000000e+00 5.792085e-16 7.888847e-34
#>  [19,] 1.000000e+00 1.838770e-12 7.077675e-28
#>  [20,] 1.000000e+00 1.007591e-15 4.932475e-33
#>  [21,] 1.000000e+00 8.003305e-13 2.988659e-29
#>  [22,] 1.000000e+00 3.239263e-15 2.212908e-32
#>  [23,] 1.000000e+00 4.293936e-21 1.455266e-41
#>  [24,] 1.000000e+00 1.684971e-12 6.043532e-29
#>  [25,] 1.000000e+00 8.089809e-13 1.237726e-29
#>  [26,] 1.000000e+00 5.879016e-14 9.233078e-32
#>  [27,] 1.000000e+00 4.281258e-14 2.980294e-31
#>  [28,] 1.000000e+00 3.431956e-15 1.078274e-32
#>  [29,] 1.000000e+00 8.148416e-16 9.369497e-34
#>  [30,] 1.000000e+00 3.725137e-15 2.067596e-33
#>  [31,] 1.000000e+00 1.082762e-14 8.120565e-33
#>  [32,] 1.000000e+00 7.954596e-14 1.340679e-30
#>  [33,] 1.000000e+00 1.170395e-16 5.229103e-34
#>  [34,] 1.000000e+00 1.545756e-16 1.980569e-33
#>  [35,] 1.000000e+00 3.086285e-15 1.505679e-33
#>  [36,] 1.000000e+00 1.316358e-17 8.441335e-37
#>  [37,] 1.000000e+00 5.208699e-16 1.073153e-33
#>  [38,] 1.000000e+00 2.498955e-17 6.344799e-36
#>  [39,] 1.000000e+00 4.915683e-18 4.121235e-38
#>  [40,] 1.000000e+00 2.858709e-15 5.053241e-33
#>  [41,] 1.000000e+00 4.731249e-17 1.745341e-35
#>  [42,] 1.000000e+00 4.189706e-16 3.337718e-36
#>  [43,] 1.000000e+00 2.030322e-18 2.239061e-38
#>  [44,] 1.000000e+00 1.174537e-13 2.402226e-30
#>  [45,] 1.000000e+00 3.840153e-12 9.692964e-28
#>  [46,] 1.000000e+00 8.104734e-16 1.494506e-34
#>  [47,] 1.000000e+00 3.196024e-15 2.328846e-32
#>  [48,] 1.000000e+00 4.643069e-17 2.929807e-36
#>  [49,] 1.000000e+00 2.647914e-15 1.695928e-32
#>  [50,] 1.000000e+00 3.633460e-16 1.516762e-34
#>  [51,] 4.414980e-22 3.552072e-01 6.447928e-01
#>  [52,] 1.058716e-18 9.643661e-01 3.563387e-02
#>  [53,] 1.636093e-24 9.729657e-02 9.027034e-01
#>  [54,] 2.934788e-13 9.999989e-01 1.061514e-06
#>  [55,] 1.418016e-20 9.269163e-01 7.308372e-02
#>  [56,] 6.343133e-17 9.996119e-01 3.881470e-04
#>  [57,] 2.880441e-20 8.010121e-01 1.989879e-01
#>  [58,] 8.884593e-05 9.999112e-01 4.919936e-11
#>  [59,] 5.157533e-20 9.481313e-01 5.186869e-02
#>  [60,] 3.557001e-11 9.999997e-01 3.315240e-07
#>  [61,] 1.889301e-07 9.999998e-01 2.501339e-10
#>  [62,] 2.902400e-15 9.997677e-01 2.323394e-04
#>  [63,] 7.313376e-14 9.999980e-01 2.017871e-06
#>  [64,] 8.929576e-20 9.771962e-01 2.280376e-02
#>  [65,] 4.113055e-09 9.999999e-01 1.162238e-07
#>  [66,] 1.420914e-18 9.695112e-01 3.048878e-02
#>  [67,] 6.717402e-17 9.991564e-01 8.435834e-04
#>  [68,] 3.551802e-13 9.999960e-01 3.972232e-06
#>  [69,] 4.586821e-20 9.961084e-01 3.891551e-03
#>  [70,] 1.064136e-11 9.999996e-01 3.545649e-07
#>  [71,] 8.183483e-21 8.135526e-01 1.864474e-01
#>  [72,] 6.299772e-14 9.999708e-01 2.918159e-05
#>  [73,] 3.889743e-23 7.726455e-01 2.273545e-01
#>  [74,] 2.487187e-19 9.921199e-01 7.880054e-03
#>  [75,] 5.333824e-17 9.983183e-01 1.681674e-03
#>  [76,] 1.728900e-18 9.826147e-01 1.738533e-02
#>  [77,] 5.532435e-23 4.900301e-01 5.099699e-01
#>  [78,] 6.253220e-26 4.638985e-02 9.536101e-01
#>  [79,] 3.533406e-18 9.957718e-01 4.228244e-03
#>  [80,] 3.378235e-08 1.000000e+00 1.222238e-08
#>  [81,] 8.372344e-11 9.999999e-01 8.373227e-08
#>  [82,] 1.133646e-09 1.000000e+00 2.195897e-08
#>  [83,] 3.573559e-12 9.999982e-01 1.792814e-06
#>  [84,] 4.604035e-24 5.142138e-01 4.857862e-01
#>  [85,] 2.345187e-16 9.996486e-01 3.514185e-04
#>  [86,] 1.545648e-17 9.866927e-01 1.330734e-02
#>  [87,] 1.499551e-21 5.950215e-01 4.049785e-01
#>  [88,] 1.071821e-18 9.988635e-01 1.136484e-03
#>  [89,] 5.293950e-13 9.999898e-01 1.023051e-05
#>  [90,] 7.105513e-13 9.999986e-01 1.396314e-06
#>  [91,] 1.238318e-15 9.999678e-01 3.223258e-05
#>  [92,] 9.214775e-19 9.888907e-01 1.110928e-02
#>  [93,] 3.504393e-13 9.999963e-01 3.724099e-06
#>  [94,] 3.056837e-05 9.999694e-01 6.648361e-11
#>  [95,] 2.144228e-14 9.999838e-01 1.615563e-05
#>  [96,] 8.934438e-14 9.999764e-01 2.359935e-05
#>  [97,] 2.779082e-14 9.999671e-01 3.293310e-05
#>  [98,] 1.863063e-16 9.992991e-01 7.008920e-04
#>  [99,] 5.368673e-03 9.946313e-01 1.595603e-11
#> [100,] 1.170520e-13 9.999879e-01 1.205306e-05
#> [101,] 9.342372e-40 7.331981e-07 9.999993e-01
#> [102,] 1.355795e-24 3.827313e-01 6.172687e-01
#> [103,] 1.415993e-40 5.197756e-07 9.999995e-01
#> [104,] 7.078729e-32 1.098448e-03 9.989016e-01
#> [105,] 3.935241e-37 1.071799e-05 9.999893e-01
#> [106,] 3.078989e-51 1.334794e-10 1.000000e+00
#> [107,] 1.370446e-16 9.999492e-01 5.076321e-05
#> [108,] 7.578449e-45 3.159087e-08 1.000000e+00
#> [109,] 1.221627e-36 5.811402e-05 9.999419e-01
#> [110,] 1.046280e-44 3.956327e-09 1.000000e+00
#> [111,] 1.775399e-27 9.010916e-03 9.909891e-01
#> [112,] 1.518738e-29 7.823942e-03 9.921761e-01
#> [113,] 2.041398e-34 6.231995e-05 9.999377e-01
#> [114,] 5.660664e-24 6.529490e-01 3.470510e-01
#> [115,] 7.166625e-27 4.895157e-02 9.510484e-01
#> [116,] 5.879819e-31 6.051577e-04 9.993948e-01
#> [117,] 1.789429e-31 9.501963e-04 9.990498e-01
#> [118,] 2.364301e-52 7.544743e-12 1.000000e+00
#> [119,] 7.542749e-57 4.302226e-12 1.000000e+00
#> [120,] 1.181638e-23 8.890101e-01 1.109899e-01
#> [121,] 4.869149e-38 2.102690e-06 9.999979e-01
#> [122,] 3.284881e-22 8.215690e-01 1.784310e-01
#> [123,] 1.223858e-52 7.611372e-11 1.000000e+00
#> [124,] 5.342028e-24 3.865481e-01 6.134519e-01
#> [125,] 6.054033e-36 1.127952e-05 9.999887e-01
#> [126,] 2.085175e-40 4.388139e-07 9.999996e-01
#> [127,] 1.763260e-22 6.698002e-01 3.301998e-01
#> [128,] 9.097872e-23 5.007794e-01 4.992206e-01
#> [129,] 4.988533e-34 1.984470e-04 9.998016e-01
#> [130,] 3.019952e-37 8.392597e-06 9.999916e-01
#> [131,] 1.422424e-43 8.290819e-08 9.999999e-01
#> [132,] 1.174192e-48 1.088190e-10 1.000000e+00
#> [133,] 1.508646e-34 1.239971e-04 9.998760e-01
#> [134,] 1.254623e-24 2.841190e-01 7.158810e-01
#> [135,] 2.771605e-29 2.548075e-02 9.745192e-01
#> [136,] 9.024027e-47 2.580492e-09 1.000000e+00
#> [137,] 2.490666e-34 3.296144e-05 9.999670e-01
#> [138,] 7.025706e-31 1.283507e-03 9.987165e-01
#> [139,] 1.751245e-21 7.873985e-01 2.126015e-01
#> [140,] 1.493664e-33 8.353069e-05 9.999165e-01
#> [141,] 1.420120e-36 8.621217e-06 9.999914e-01
#> [142,] 5.196881e-31 4.407646e-04 9.995592e-01
#> [143,] 1.355795e-24 3.827313e-01 6.172687e-01
#> [144,] 5.782379e-40 5.741210e-07 9.999994e-01
#> [145,] 5.062652e-38 1.718831e-06 9.999983e-01
#> [146,] 2.055977e-31 5.095741e-04 9.994904e-01
#> [147,] 7.527149e-26 1.785571e-01 8.214429e-01
#> [148,] 6.202938e-29 4.995365e-03 9.950046e-01
#> [149,] 5.809956e-31 4.638463e-04 9.995362e-01
#> [150,] 4.398339e-24 2.979812e-01 7.020188e-01
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
#> [1] "EII"
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
#> Sepal.Length     0.148829    0.000000     0.000000    0.000000
#> Sepal.Width      0.000000    0.148829     0.000000    0.000000
#> Petal.Length     0.000000    0.000000     0.148829    0.000000
#> Petal.Width      0.000000    0.000000     0.000000    0.148829
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.148829    0.000000     0.000000    0.000000
#> Sepal.Width      0.000000    0.148829     0.000000    0.000000
#> Petal.Length     0.000000    0.000000     0.148829    0.000000
#> Petal.Width      0.000000    0.000000     0.000000    0.148829
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.148829    0.000000     0.000000    0.000000
#> Sepal.Width      0.000000    0.148829     0.000000    0.000000
#> Petal.Length     0.000000    0.000000     0.148829    0.000000
#> Petal.Width      0.000000    0.000000     0.000000    0.148829
#> 
#> 
#> $parameters$variance$Sigma
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.148829    0.000000     0.000000    0.000000
#> Sepal.Width      0.000000    0.148829     0.000000    0.000000
#> Petal.Length     0.000000    0.000000     0.148829    0.000000
#> Petal.Width      0.000000    0.000000     0.000000    0.148829
#> 
#> $parameters$variance$sigmasq
#> [1] 0.148829
#> 
#> $parameters$variance$scale
#> [1] 0.148829
#> 
#> 
#> 
#> $loglik
#> [1] -414.698
#> 
#> attr(,"returnCode")
#> [1] 0
```
