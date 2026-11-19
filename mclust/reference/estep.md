# E-step for parameterized Gaussian mixture models.

Implements the expectation step of EM algorithm for parameterized
Gaussian mixture models.

## Usage

``` r
estep(data, modelName, parameters, warn = NULL, ...)
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

- parameters:

  A names list giving the parameters of the model. The components are as
  follows:

  `pro`

  : Mixing proportions for the components of the mixture. If the model
    includes a Poisson term for noise, there should be one more mixing
    proportion than the number of Gaussian components.

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

  : An estimate of the reciprocal hypervolume of the data region. If set
    to NULL or a negative value, the default is determined by applying
    function `hypvol` to the data. Used only when `pro` includes an
    additional mixing proportion for a noise component.

- warn:

  A logical value indicating whether or not a warning should be issued
  when computations fail. The default is `warn=FALSE`.

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

  The input parameters.

- loglik:

  The log-likelihood for the data in the mixture model.

- Attributes:

  `"WARNING"`: an appropriate warning if problems are encountered in the
  computations.

## See also

[`estepE`](https://mclust-org.github.io/mclust/reference/estepE.md),
...,
[`estepVVV`](https://mclust-org.github.io/mclust/reference/estepE.md),
[`em`](https://mclust-org.github.io/mclust/reference/em.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`mclust.options`](https://mclust-org.github.io/mclust/reference/mclust.options.md)
[`mclustVariance`](https://mclust-org.github.io/mclust/reference/mclustVariance.md)

## Examples

``` r
# \donttest{
msEst <- mstep(modelName = "VVV", data = iris[,-5], z = unmap(iris[,5]))
names(msEst)
#> [1] "modelName"  "prior"      "n"          "d"          "G"         
#> [6] "z"          "parameters"

estep(modelName = msEst$modelName, data = iris[,-5],
      parameters = msEst$parameters)# }
#> $modelName
#> [1] "VVV"
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
#>   [1,]  1.000000e+00 1.531298e-26 4.631660e-42
#>   [2,]  1.000000e+00 3.341913e-19 2.782813e-35
#>   [3,]  1.000000e+00 5.970550e-22 6.659211e-37
#>   [4,]  1.000000e+00 3.629271e-19 2.069429e-32
#>   [5,]  1.000000e+00 9.920165e-28 3.097629e-42
#>   [6,]  1.000000e+00 4.473024e-27 2.035241e-41
#>   [7,]  1.000000e+00 1.023721e-21 5.400089e-35
#>   [8,]  1.000000e+00 3.951048e-24 1.188427e-38
#>   [9,]  1.000000e+00 6.348304e-17 2.222620e-30
#>  [10,]  1.000000e+00 3.282086e-21 1.202430e-35
#>  [11,]  1.000000e+00 1.619962e-29 2.944164e-45
#>  [12,]  1.000000e+00 7.891556e-23 2.286635e-35
#>  [13,]  1.000000e+00 2.172110e-20 1.922311e-35
#>  [14,]  1.000000e+00 1.830223e-20 8.325932e-35
#>  [15,]  1.000000e+00 1.347076e-37 5.789535e-57
#>  [16,]  1.000000e+00 2.707999e-38 1.150952e-53
#>  [17,]  1.000000e+00 1.194101e-30 1.742225e-47
#>  [18,]  1.000000e+00 6.948699e-25 1.025185e-40
#>  [19,]  1.000000e+00 5.504880e-28 4.538366e-44
#>  [20,]  1.000000e+00 2.069605e-28 1.972755e-42
#>  [21,]  1.000000e+00 7.479277e-23 3.988370e-38
#>  [22,]  1.000000e+00 5.032153e-25 1.181412e-39
#>  [23,]  1.000000e+00 7.912651e-28 1.053516e-42
#>  [24,]  1.000000e+00 5.071504e-16 1.408708e-30
#>  [25,]  1.000000e+00 5.079161e-20 5.211806e-30
#>  [26,]  1.000000e+00 9.462960e-18 6.946195e-33
#>  [27,]  1.000000e+00 8.007503e-20 3.118746e-34
#>  [28,]  1.000000e+00 5.238689e-26 2.393533e-41
#>  [29,]  1.000000e+00 2.981616e-25 9.025368e-42
#>  [30,]  1.000000e+00 6.305750e-20 1.248894e-32
#>  [31,]  1.000000e+00 7.474435e-19 1.745401e-32
#>  [32,]  1.000000e+00 1.143232e-21 6.109943e-39
#>  [33,]  1.000000e+00 1.712290e-37 3.848459e-49
#>  [34,]  1.000000e+00 4.783083e-39 6.903779e-54
#>  [35,]  1.000000e+00 8.338572e-20 8.305829e-35
#>  [36,]  1.000000e+00 2.989188e-23 1.006538e-40
#>  [37,]  1.000000e+00 6.887941e-28 1.135265e-46
#>  [38,]  1.000000e+00 3.867251e-29 2.455633e-42
#>  [39,]  1.000000e+00 1.417061e-18 2.036480e-32
#>  [40,]  1.000000e+00 2.473260e-24 1.859072e-39
#>  [41,]  1.000000e+00 3.020988e-25 3.902559e-41
#>  [42,]  1.000000e+00 2.560802e-10 2.997338e-26
#>  [43,]  1.000000e+00 7.792261e-21 5.097178e-34
#>  [44,]  1.000000e+00 4.049190e-17 8.149856e-31
#>  [45,]  1.000000e+00 2.350953e-22 5.235460e-34
#>  [46,]  1.000000e+00 1.246886e-17 3.033197e-33
#>  [47,]  1.000000e+00 1.560647e-29 3.679695e-42
#>  [48,]  1.000000e+00 4.767026e-21 1.084422e-34
#>  [49,]  1.000000e+00 2.532098e-29 1.912018e-44
#>  [50,]  1.000000e+00 1.942649e-23 4.162666e-39
#>  [51,]  4.427741e-92 9.999635e-01 3.651562e-05
#>  [52,]  1.492993e-83 9.996963e-01 3.036808e-04
#>  [53,] 2.323794e-104 9.985558e-01 1.444189e-03
#>  [54,]  3.820558e-64 9.974857e-01 2.514348e-03
#>  [55,]  2.529967e-92 9.975986e-01 2.401407e-03
#>  [56,]  3.438709e-79 9.895919e-01 1.040812e-02
#>  [57,]  9.128977e-94 9.951481e-01 4.851888e-03
#>  [58,]  5.148001e-34 9.999951e-01 4.872649e-06
#>  [59,]  1.880560e-86 9.998169e-01 1.830687e-04
#>  [60,]  1.542139e-60 9.941598e-01 5.840159e-03
#>  [61,]  2.992578e-42 9.999480e-01 5.203024e-05
#>  [62,]  1.708982e-72 9.987050e-01 1.295008e-03
#>  [63,]  3.769219e-61 9.999887e-01 1.133958e-05
#>  [64,]  1.085898e-90 9.892610e-01 1.073895e-02
#>  [65,]  9.208897e-47 9.999914e-01 8.648146e-06
#>  [66,]  4.247145e-79 9.999894e-01 1.063020e-05
#>  [67,]  1.205637e-83 9.747235e-01 2.527649e-02
#>  [68,]  6.402952e-58 9.998866e-01 1.133547e-04
#>  [69,]  5.514687e-92 8.146259e-01 1.853741e-01
#>  [70,]  9.644338e-55 9.999704e-01 2.957058e-05
#>  [71,] 8.144832e-106 3.284513e-01 6.715487e-01
#>  [72,]  5.073314e-62 9.999918e-01 8.225435e-06
#>  [73,] 4.134551e-107 6.987624e-01 3.012376e-01
#>  [74,]  4.817238e-86 9.735004e-01 2.649958e-02
#>  [75,]  2.448343e-73 9.999831e-01 1.685110e-05
#>  [76,]  9.080831e-80 9.999734e-01 2.661419e-05
#>  [77,] 3.308249e-100 9.986285e-01 1.371485e-03
#>  [78,] 6.162406e-115 8.630616e-01 1.369384e-01
#>  [79,]  3.606949e-85 9.927366e-01 7.263396e-03
#>  [80,]  1.827587e-39 9.999999e-01 6.520517e-08
#>  [81,]  5.062737e-52 9.999754e-01 2.463634e-05
#>  [82,]  1.118752e-46 9.999953e-01 4.706239e-06
#>  [83,]  5.410273e-56 9.999908e-01 9.231469e-06
#>  [84,] 1.930587e-116 1.473576e-01 8.526424e-01
#>  [85,]  1.417725e-83 9.454001e-01 5.459993e-02
#>  [86,]  4.017430e-84 9.963529e-01 3.647118e-03
#>  [87,]  2.333373e-94 9.994678e-01 5.322498e-04
#>  [88,]  4.416063e-83 9.990213e-01 9.787027e-04
#>  [89,]  1.824762e-62 9.997972e-01 2.027893e-04
#>  [90,]  9.617210e-63 9.990635e-01 9.364656e-04
#>  [91,]  1.267274e-73 9.817631e-01 1.823689e-02
#>  [92,]  1.366483e-85 9.972534e-01 2.746587e-03
#>  [93,]  2.483385e-60 9.999664e-01 3.362787e-05
#>  [94,]  8.736184e-35 9.999970e-01 2.980291e-06
#>  [95,]  1.632088e-68 9.987856e-01 1.214423e-03
#>  [96,]  3.007660e-63 9.997647e-01 2.352778e-04
#>  [97,]  4.649893e-67 9.996298e-01 3.702063e-04
#>  [98,]  2.213572e-72 9.999435e-01 5.649197e-05
#>  [99,]  4.667070e-28 9.999995e-01 5.312475e-07
#> [100,]  2.415298e-64 9.997932e-01 2.068042e-04
#> [101,] 5.431127e-203 2.210439e-09 1.000000e+00
#> [102,] 1.642385e-128 3.822927e-04 9.996177e-01
#> [103,] 2.935239e-181 4.340914e-05 9.999566e-01
#> [104,] 1.281149e-148 5.160913e-03 9.948391e-01
#> [105,] 8.739248e-178 1.967500e-06 9.999980e-01
#> [106,] 1.169524e-228 1.454416e-06 9.999985e-01
#> [107,]  3.059959e-95 3.396516e-03 9.966035e-01
#> [108,] 1.199215e-196 3.541718e-05 9.999646e-01
#> [109,] 1.063382e-166 1.059411e-04 9.998941e-01
#> [110,] 1.856330e-207 4.601010e-08 1.000000e+00
#> [111,] 1.697668e-129 5.349751e-03 9.946502e-01
#> [112,] 8.800846e-140 8.088253e-04 9.991912e-01
#> [113,] 9.146951e-158 5.195774e-05 9.999480e-01
#> [114,] 1.691089e-130 1.177827e-06 9.999988e-01
#> [115,] 1.529536e-156 1.686232e-13 1.000000e+00
#> [116,] 8.247834e-156 4.558397e-08 1.000000e+00
#> [117,] 6.788928e-143 3.037798e-02 9.696220e-01
#> [118,] 2.930198e-228 2.149645e-04 9.997850e-01
#> [119,] 9.079687e-265 8.879452e-11 1.000000e+00
#> [120,] 2.303150e-113 3.790989e-02 9.620901e-01
#> [121,] 3.843003e-177 3.783066e-07 9.999996e-01
#> [122,] 9.595809e-124 1.976052e-05 9.999802e-01
#> [123,] 9.581470e-236 2.100380e-07 9.999998e-01
#> [124,] 1.065022e-115 2.567755e-02 9.743225e-01
#> [125,] 1.850244e-164 7.609181e-04 9.992391e-01
#> [126,] 1.221324e-172 6.521046e-03 9.934790e-01
#> [127,] 7.743971e-110 5.265561e-02 9.473444e-01
#> [128,] 3.432568e-112 1.442112e-01 8.557888e-01
#> [129,] 1.151967e-163 4.572452e-06 9.999954e-01
#> [130,] 3.706143e-157 1.797263e-02 9.820274e-01
#> [131,] 6.470511e-190 1.243961e-04 9.998756e-01
#> [132,] 3.730831e-201 8.236967e-03 9.917630e-01
#> [133,] 6.157428e-169 1.310745e-07 9.999999e-01
#> [134,] 2.506178e-113 6.022880e-01 3.977120e-01
#> [135,] 3.266558e-138 1.780195e-04 9.998220e-01
#> [136,] 2.569119e-207 8.589823e-08 9.999999e-01
#> [137,] 8.687854e-175 3.751507e-08 1.000000e+00
#> [138,] 2.183741e-141 4.651619e-02 9.534838e-01
#> [139,] 1.312185e-107 1.339905e-01 8.660095e-01
#> [140,] 9.141262e-152 1.555449e-04 9.998445e-01
#> [141,] 2.098453e-178 6.704480e-10 1.000000e+00
#> [142,] 2.368114e-148 4.059502e-09 1.000000e+00
#> [143,] 1.642385e-128 3.822927e-04 9.996177e-01
#> [144,] 1.224376e-187 6.897983e-07 9.999993e-01
#> [145,] 8.404700e-188 1.144473e-10 1.000000e+00
#> [146,] 1.099698e-153 2.083183e-09 1.000000e+00
#> [147,] 5.133200e-127 1.360843e-04 9.998639e-01
#> [148,] 1.472988e-136 9.106700e-04 9.990893e-01
#> [149,] 1.208779e-158 8.694131e-07 9.999991e-01
#> [150,] 2.673436e-121 5.663609e-02 9.433639e-01
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
#> Sepal.Length     0.261104     0.08348      0.17924    0.054664
#> Sepal.Width      0.083480     0.09650      0.08100    0.040380
#> Petal.Length     0.179240     0.08100      0.21640    0.071640
#> Petal.Width      0.054664     0.04038      0.07164    0.038324
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length     0.396256    0.091888     0.297224    0.048112
#> Sepal.Width      0.091888    0.101924     0.069952    0.046676
#> Petal.Length     0.297224    0.069952     0.298496    0.047848
#> Petal.Width      0.048112    0.046676     0.047848    0.073924
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
#> Sepal.Length    0.5109834   0.1633713   0.35077463  0.10697804
#> Sepal.Width     0.0000000  -0.2642155  -0.08967492 -0.08668251
#> Petal.Length    0.0000000   0.0000000   0.29208829  0.09018359
#> Petal.Width     0.0000000   0.0000000   0.00000000 -0.10598473
#> 
#> , , 3
#> 
#>              Sepal.Length Sepal.Width Petal.Length Petal.Width
#> Sepal.Length   -0.6294887  -0.1459724 -0.472167346 -0.07643029
#> Sepal.Width     0.0000000  -0.2839297 -0.003622656 -0.12509889
#> Petal.Length    0.0000000   0.0000000  0.274847001  0.04113898
#> Petal.Width     0.0000000   0.0000000  0.000000000 -0.22525600
#> 
#> 
#> 
#> 
#> $loglik
#> [1] -182.9208
#> 
#> attr(,"returnCode")
#> [1] 0
```
