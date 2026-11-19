# EM algorithm starting with M-step for a parameterized Gaussian mixture model

Implements the EM algorithm for a parameterized Gaussian mixture model,
starting with the maximization step.

## Usage

``` r
meE(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meV(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meX(data, prior = NULL, warn = NULL, ...)
meEII(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVII(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meEEI(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVEI(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meEVI(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVVI(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meEEE(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVEE(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meEVE(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVVE(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meEEV(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVEV(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meEVV(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meVVV(data, z, prior=NULL, control=emControl(), Vinv=NULL, warn=NULL, ...)
meXII(data, prior = NULL, warn = NULL, ...)
meXXI(data, prior = NULL, warn = NULL, ...)
meXXX(data, prior = NULL, warn = NULL, ...)
```

## Arguments

- data:

  A numeric vector, matrix, or data frame of observations. Categorical
  variables are not allowed. If a matrix or data frame, rows correspond
  to observations and columns correspond to variables.

- z:

  A matrix whose `[i,k]`th entry is the conditional probability of the
  ith observation belonging to the *k*th component of the mixture.

- prior:

  Specification of a conjugate prior on the means and variances. The
  default assumes no prior.

- control:

  A list of control parameters for EM. The defaults are set by the call
  [`emControl()`](https://mclust-org.github.io/mclust/reference/emControl.md).

- Vinv:

  An estimate of the reciprocal hypervolume of the data region, when the
  model is to include a noise term. Set to a negative value or zero if a
  noise term is desired, but an estimate is unavailable — in that case
  function `hypvol` will be used to obtain the estimate. The default is
  not to assume a noise term in the model through the setting
  `Vinv=NULL`.

- warn:

  A logical value indicating whether or not certain warnings (usually
  related to singularity) should be issued when the estimation fails.
  The default is given by `mclust.options("warn")`.

- ... :

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
  computations.\

## See also

[`em`](https://mclust-org.github.io/mclust/reference/em.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md)

## Examples

``` r
meVVV(data = iris[,-5], z = unmap(iris[,5]))
#> $modelName
#> [1] "VVV"
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
#>                 [,1]         [,2]         [,3]
#>   [1,]  1.000000e+00 1.340380e-44 1.861339e-34
#>   [2,]  1.000000e+00 2.201405e-31 6.676298e-28
#>   [3,]  1.000000e+00 1.896748e-36 1.102178e-29
#>   [4,]  1.000000e+00 3.488647e-32 6.409600e-26
#>   [5,]  1.000000e+00 4.393475e-47 7.745885e-35
#>   [6,]  1.000000e+00 1.278514e-45 9.141846e-35
#>   [7,]  1.000000e+00 1.725033e-36 1.528128e-28
#>   [8,]  1.000000e+00 1.013323e-40 1.687173e-31
#>   [9,]  1.000000e+00 6.118503e-28 6.204452e-24
#>  [10,]  1.000000e+00 3.941874e-36 2.494386e-28
#>  [11,]  1.000000e+00 4.448705e-50 1.929162e-37
#>  [12,]  1.000000e+00 3.292550e-39 7.578657e-29
#>  [13,]  1.000000e+00 2.565267e-34 5.452124e-28
#>  [14,]  1.000000e+00 9.106612e-35 1.909555e-27
#>  [15,]  1.000000e+00 3.125165e-63 1.593148e-47
#>  [16,]  1.000000e+00 1.896077e-64 4.876149e-46
#>  [17,]  1.000000e+00 5.596191e-50 1.021588e-39
#>  [18,]  1.000000e+00 6.308031e-41 2.149597e-33
#>  [19,]  1.000000e+00 2.370660e-47 1.329559e-36
#>  [20,]  1.000000e+00 4.052147e-48 1.898553e-35
#>  [21,]  1.000000e+00 3.784797e-39 8.212299e-31
#>  [22,]  1.000000e+00 2.257628e-41 6.202577e-33
#>  [23,]  1.000000e+00 6.419874e-48 7.020340e-35
#>  [24,]  1.000000e+00 6.639759e-26 1.443374e-24
#>  [25,]  1.000000e+00 7.175865e-36 2.846195e-24
#>  [26,]  1.000000e+00 1.227836e-29 6.692266e-26
#>  [27,]  1.000000e+00 1.658055e-32 6.789689e-28
#>  [28,]  1.000000e+00 6.872462e-44 7.516822e-34
#>  [29,]  1.000000e+00 4.516919e-42 6.057310e-34
#>  [30,]  1.000000e+00 6.438483e-34 2.970895e-26
#>  [31,]  1.000000e+00 9.491263e-32 6.566490e-26
#>  [32,]  1.000000e+00 1.359677e-34 1.377621e-31
#>  [33,]  1.000000e+00 1.807081e-66 1.923514e-41
#>  [34,]  1.000000e+00 2.640245e-67 1.082266e-45
#>  [35,]  1.000000e+00 6.866910e-33 9.278292e-28
#>  [36,]  1.000000e+00 5.302053e-38 1.312420e-32
#>  [37,]  1.000000e+00 3.083368e-46 5.787644e-38
#>  [38,]  1.000000e+00 1.102357e-50 7.771785e-35
#>  [39,]  1.000000e+00 8.786404e-31 1.026566e-25
#>  [40,]  1.000000e+00 6.077402e-41 4.042171e-32
#>  [41,]  1.000000e+00 1.877063e-41 1.037298e-33
#>  [42,]  1.000000e+00 3.904478e-15 2.495765e-19
#>  [43,]  1.000000e+00 2.822962e-35 2.367557e-27
#>  [44,]  1.000000e+00 1.608040e-27 5.055302e-25
#>  [45,]  1.000000e+00 3.082923e-39 2.560673e-28
#>  [46,]  1.000000e+00 4.000241e-28 2.471836e-26
#>  [47,]  1.000000e+00 2.987492e-51 3.564570e-35
#>  [48,]  1.000000e+00 2.344766e-35 6.084895e-28
#>  [49,]  1.000000e+00 7.194407e-50 8.191014e-37
#>  [50,]  1.000000e+00 5.327318e-39 1.230925e-31
#>  [51,]  3.032539e-92 9.997227e-01 2.773164e-04
#>  [52,]  1.150038e-83 9.986288e-01 1.371182e-03
#>  [53,] 1.797174e-104 9.944941e-01 5.505859e-03
#>  [54,]  1.115012e-63 9.320826e-01 6.791742e-02
#>  [55,]  4.881945e-92 9.705411e-01 2.945894e-02
#>  [56,]  2.954459e-79 9.685784e-01 3.142162e-02
#>  [57,]  8.550867e-94 9.865626e-01 1.343744e-02
#>  [58,]  3.765759e-34 9.998441e-01 1.558591e-04
#>  [59,]  1.336478e-86 9.986175e-01 1.382503e-03
#>  [60,]  1.747387e-60 9.669572e-01 3.304278e-02
#>  [61,]  3.021885e-42 9.982058e-01 1.794178e-03
#>  [62,]  1.635100e-72 9.927992e-01 7.200848e-03
#>  [63,]  4.079064e-61 9.991369e-01 8.630789e-04
#>  [64,]  1.012713e-90 9.670966e-01 3.290342e-02
#>  [65,]  7.559543e-47 9.997995e-01 2.005120e-04
#>  [66,]  3.050597e-79 9.998767e-01 1.232632e-04
#>  [67,]  1.489175e-83 9.246757e-01 7.532431e-02
#>  [68,]  8.933659e-58 9.965960e-01 3.404046e-03
#>  [69,]  3.934648e-92 2.925580e-03 9.970744e-01
#>  [70,]  7.242827e-55 9.995024e-01 4.975685e-04
#>  [71,] 7.707744e-106 5.397566e-02 9.460243e-01
#>  [72,]  3.850909e-62 9.998443e-01 1.556841e-04
#>  [73,] 4.713812e-107 4.325311e-02 9.567469e-01
#>  [74,]  3.438056e-86 9.136402e-01 8.635976e-02
#>  [75,]  1.817276e-73 9.997824e-01 2.176209e-04
#>  [76,]  6.615926e-80 9.997244e-01 2.755883e-04
#>  [77,] 3.767869e-100 9.879246e-01 1.207544e-02
#>  [78,] 1.254096e-114 3.361880e-01 6.638120e-01
#>  [79,]  4.845433e-85 9.646654e-01 3.533460e-02
#>  [80,]  2.507098e-39 9.999864e-01 1.363479e-05
#>  [81,]  3.878089e-52 9.995100e-01 4.900484e-04
#>  [82,]  8.943244e-47 9.998196e-01 1.803552e-04
#>  [83,]  4.112093e-56 9.998282e-01 1.717912e-04
#>  [84,] 1.761949e-116 6.973648e-03 9.930264e-01
#>  [85,]  1.873570e-83 8.494274e-01 1.505726e-01
#>  [86,]  4.330051e-84 9.873428e-01 1.265717e-02
#>  [87,]  1.827651e-94 9.975530e-01 2.447044e-03
#>  [88,]  2.149510e-82 9.204475e-01 7.955247e-02
#>  [89,]  1.777598e-62 9.985437e-01 1.456316e-03
#>  [90,]  1.223594e-62 9.899298e-01 1.007019e-02
#>  [91,]  1.057937e-73 9.412283e-01 5.877173e-02
#>  [92,]  1.119977e-85 9.910150e-01 8.985033e-03
#>  [93,]  1.917546e-60 9.995110e-01 4.889726e-04
#>  [94,]  5.809300e-35 9.998832e-01 1.168043e-04
#>  [95,]  1.439348e-68 9.937569e-01 6.243104e-03
#>  [96,]  4.058872e-63 9.973306e-01 2.669419e-03
#>  [97,]  3.805944e-67 9.979620e-01 2.037953e-03
#>  [98,]  1.658725e-72 9.994971e-01 5.028802e-04
#>  [99,]  3.308948e-28 9.999418e-01 5.818736e-05
#> [100,]  1.884676e-64 9.986048e-01 1.395245e-03
#> [101,] 3.010949e-203 5.530199e-17 1.000000e+00
#> [102,] 1.753436e-128 1.077534e-07 9.999999e-01
#> [103,] 3.185553e-181 2.653929e-09 1.000000e+00
#> [104,] 1.348868e-148 4.312002e-05 9.999569e-01
#> [105,] 9.131952e-178 7.830553e-12 1.000000e+00
#> [106,] 1.376312e-228 2.727688e-11 1.000000e+00
#> [107,]  3.275601e-95 1.168041e-05 9.999883e-01
#> [108,] 1.485822e-196 1.996513e-07 9.999998e-01
#> [109,] 1.096470e-166 5.615666e-09 1.000000e+00
#> [110,] 1.838986e-207 5.153929e-13 1.000000e+00
#> [111,] 1.407473e-129 6.017637e-05 9.999398e-01
#> [112,] 7.916908e-140 1.650636e-07 9.999998e-01
#> [113,] 9.136916e-158 5.957065e-09 1.000000e+00
#> [114,] 1.956477e-130 2.233857e-12 1.000000e+00
#> [115,] 2.022253e-156 3.404704e-22 1.000000e+00
#> [116,] 9.392240e-156 1.026315e-12 1.000000e+00
#> [117,] 6.602144e-143 7.801735e-04 9.992198e-01
#> [118,] 2.748870e-228 5.421043e-06 9.999946e-01
#> [119,] 7.429637e-265 3.265440e-22 1.000000e+00
#> [120,] 1.370738e-113 5.282640e-05 9.999472e-01
#> [121,] 4.351788e-177 4.577282e-12 1.000000e+00
#> [122,] 1.054301e-123 3.279642e-09 1.000000e+00
#> [123,] 1.162139e-235 1.885723e-12 1.000000e+00
#> [124,] 6.201270e-116 5.216528e-05 9.999478e-01
#> [125,] 1.891890e-164 2.001713e-06 9.999980e-01
#> [126,] 1.323823e-172 7.420387e-04 9.992580e-01
#> [127,] 4.617100e-110 3.526679e-04 9.996473e-01
#> [128,] 2.794272e-112 5.946968e-03 9.940530e-01
#> [129,] 1.280855e-163 1.922717e-11 1.000000e+00
#> [130,] 3.061896e-157 4.726577e-03 9.952734e-01
#> [131,] 6.533587e-190 3.218376e-08 1.000000e+00
#> [132,] 4.771065e-201 6.006480e-03 9.939935e-01
#> [133,] 6.908783e-169 3.468057e-14 1.000000e+00
#> [134,] 2.744809e-113 2.199157e-01 7.800843e-01
#> [135,] 3.663047e-138 3.079362e-05 9.999692e-01
#> [136,] 2.609736e-207 1.978334e-14 1.000000e+00
#> [137,] 7.627542e-175 3.636807e-13 1.000000e+00
#> [138,] 2.136992e-141 2.459565e-03 9.975404e-01
#> [139,] 1.018504e-107 5.760480e-03 9.942395e-01
#> [140,] 8.410204e-152 7.863015e-08 9.999999e-01
#> [141,] 2.436458e-178 1.184065e-16 1.000000e+00
#> [142,] 2.086983e-148 3.695558e-14 1.000000e+00
#> [143,] 1.753436e-128 1.077534e-07 9.999999e-01
#> [144,] 1.276099e-187 4.205111e-12 1.000000e+00
#> [145,] 9.196112e-188 2.313054e-17 1.000000e+00
#> [146,] 1.126391e-153 3.179633e-15 1.000000e+00
#> [147,] 3.395534e-127 2.310806e-09 1.000000e+00
#> [148,] 1.300448e-136 1.192490e-06 9.999988e-01
#> [149,] 1.141901e-158 1.374617e-10 1.000000e+00
#> [150,] 2.488135e-121 1.722216e-03 9.982778e-01
#> 
#> $parameters
#> $parameters$pro
#> [1] 0.3333333 0.2995864 0.3670803
#> 
#> $parameters$mean
#>               [,1]     [,2]     [,3]
#> Sepal.Length 5.006 5.915306 6.544949
#> Sepal.Width  3.428 2.777875 2.948818
#> Petal.Length 1.462 4.202227 5.480372
#> Petal.Width  0.246 1.297229 1.985127
#> 
#> $parameters$variance
#> $parameters$variance$modelName
#> [1] "VVV"
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
#> Sepal.Length     0.121764    0.097232     0.016028    0.010124
#> Sepal.Width      0.097232    0.140816     0.011464    0.009112
#> Petal.Length     0.016028    0.011464     0.029556    0.005948
#> Petal.Width      0.010124    0.009112     0.005948    0.010884
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.27533588  0.09687192   0.18477153  0.05444692
#> Sepal.Width    0.09687192  0.09262785   0.09111654  0.04299505
#> Petal.Length   0.18477153  0.09111654   0.20089896  0.06109388
#> Petal.Width    0.05444692  0.04299505   0.06109388  0.03205203
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   0.38707068  0.09220736   0.30268373  0.06147257
#> Sepal.Width    0.09220736  0.11034902   0.08419484  0.05595573
#> Petal.Length   0.30268373  0.08419484   0.32736949  0.07416052
#> Petal.Width    0.06147257  0.05595573   0.07416052  0.08559833
#> 
#> 
#> $parameters$variance$cholsigma
#> , , 1
#> 
#>              Sepal.Length Sepal.Width Petal.Length  Petal.Width
#> Sepal.Length    -0.348947  -0.2786440 -0.045932479 -0.029013003
#> Sepal.Width      0.000000   0.2513434 -0.005310709  0.004088826
#> Petal.Length     0.000000   0.0000000 -0.165583827 -0.028004398
#> Petal.Width      0.000000   0.0000000  0.000000000  0.096131581
#> 
#> , , 2
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length    0.5247246   0.1846148    0.3521305  0.10376286
#> Sepal.Width     0.0000000  -0.2419612   -0.1079018 -0.09852359
#> Petal.Length    0.0000000   0.0000000    0.2554609  0.05450909
#> Petal.Width     0.0000000   0.0000000    0.0000000 -0.09277479
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     -0.62215  -0.1482076  -0.48651244 -0.09880666
#> Sepal.Width       0.00000  -0.2972937  -0.04066688 -0.13895968
#> Petal.Length      0.00000   0.0000000  -0.29836444 -0.06850276
#> Petal.Width       0.00000   0.0000000   0.00000000 -0.22766895
#> 
#> 
#> 
#> $parameters$Vinv
#> NULL
#> 
#> 
#> $control
#> $control$eps
#> [1] 2.220446e-16
#> 
#> $control$tol
#> [1] 1.000000e-05 1.490116e-08
#> 
#> $control$itmax
#> [1] 2147483647 2147483647
#> 
#> $control$equalPro
#> [1] FALSE
#> 
#> 
#> $loglik
#> [1] -180.1859
#> 
#> attr(,"info")
#>   iterations        error 
#> 1.100000e+01 4.438788e-06 
#> attr(,"returnCode")
#> [1] 0
```
