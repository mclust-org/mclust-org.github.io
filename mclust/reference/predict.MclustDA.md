# Classify multivariate observations by Gaussian finite mixture modeling

Classify multivariate observations based on Gaussian finite mixture
models estimated by
[`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

## Usage

``` r
# S3 method for class 'MclustDA'
predict(object, newdata, prop = object$prop, ...)
```

## Arguments

- object:

  an object of class `'MclustDA'` resulting from a call to
  [`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

- newdata:

  a data frame or matrix giving the data. If missing the train data
  obtained from the call to
  [`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md)
  are classified.

- prop:

  the class proportions or prior class probabilities to belong to each
  class; by default, this is set at the class proportions in the
  training data.

- ...:

  further arguments passed to or from other methods.

## Value

Returns a list of with the following components:

- classification:

  a factor of predicted class labels for `newdata`.

- z:

  a matrix whose *\[i,k\]*th entry is the probability that observation
  *i* in `newdata` belongs to the *k*th class.

## Author

Luca Scrucca

## See also

[`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

## Examples

``` r
# \donttest{
odd <- seq(from = 1, to = nrow(iris), by = 2)
even <- odd + 1
X.train <- iris[odd,-5]
Class.train <- iris[odd,5]
X.test <- iris[even,-5]
Class.test <- iris[even,5]

irisMclustDA <- MclustDA(X.train, Class.train)

predTrain <- predict(irisMclustDA)
predTrain
#> $classification
#>  [1] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [7] setosa     setosa     setosa     setosa     setosa     setosa    
#> [13] setosa     setosa     setosa     setosa     setosa     setosa    
#> [19] setosa     setosa     setosa     setosa     setosa     setosa    
#> [25] setosa     versicolor versicolor versicolor versicolor versicolor
#> [31] versicolor versicolor versicolor versicolor versicolor versicolor
#> [37] versicolor versicolor versicolor versicolor versicolor versicolor
#> [43] versicolor versicolor versicolor versicolor versicolor versicolor
#> [49] versicolor versicolor virginica  virginica  virginica  virginica 
#> [55] virginica  virginica  virginica  virginica  virginica  virginica 
#> [61] virginica  virginica  virginica  virginica  virginica  virginica 
#> [67] virginica  virginica  virginica  virginica  virginica  virginica 
#> [73] virginica  virginica  virginica 
#> Levels: setosa versicolor virginica
#> 
#> $z
#>              setosa   versicolor    virginica
#>  [1,]  1.000000e+00 8.738718e-25 2.855764e-40
#>  [2,]  1.000000e+00 3.763843e-20 4.191170e-35
#>  [3,]  1.000000e+00 1.192318e-25 1.628578e-40
#>  [4,]  1.000000e+00 1.750218e-18 2.283705e-32
#>  [5,]  1.000000e+00 1.977277e-15 1.394389e-28
#>  [6,]  1.000000e+00 5.467422e-29 1.537812e-44
#>  [7,]  1.000000e+00 2.759636e-19 1.532196e-33
#>  [8,]  1.000000e+00 1.403399e-37 3.948656e-57
#>  [9,]  1.000000e+00 1.272672e-27 5.020345e-46
#> [10,]  1.000000e+00 2.998149e-28 5.201993e-44
#> [11,]  1.000000e+00 1.777956e-23 1.196557e-37
#> [12,]  1.000000e+00 2.155743e-24 1.313946e-40
#> [13,]  1.000000e+00 6.685010e-20 1.573532e-29
#> [14,]  1.000000e+00 5.386284e-17 8.640940e-32
#> [15,]  1.000000e+00 6.242576e-24 4.500542e-40
#> [16,]  1.000000e+00 7.787919e-18 7.087150e-31
#> [17,]  1.000000e+00 4.607979e-36 4.867880e-49
#> [18,]  1.000000e+00 8.821491e-19 3.074362e-33
#> [19,]  1.000000e+00 6.216587e-28 3.741545e-46
#> [20,]  1.000000e+00 5.036637e-17 9.059046e-31
#> [21,]  1.000000e+00 1.521644e-22 4.556289e-39
#> [22,]  1.000000e+00 6.593237e-19 3.086722e-32
#> [23,]  1.000000e+00 1.357970e-20 2.658524e-33
#> [24,]  1.000000e+00 2.718795e-28 1.554929e-41
#> [25,]  1.000000e+00 2.177485e-28 1.601509e-43
#> [26,]  3.958805e-98 9.999991e-01 8.952163e-07
#> [27,] 1.848038e-113 9.999889e-01 1.114528e-05
#> [28,] 5.757618e-103 9.969487e-01 3.051293e-03
#> [29,] 4.092148e-111 9.974765e-01 2.523452e-03
#> [30,]  4.075046e-88 9.999924e-01 7.631479e-06
#> [31,]  3.354121e-49 9.999910e-01 8.973135e-06
#> [32,]  1.699797e-57 9.999925e-01 7.467309e-06
#> [33,]  3.171219e-60 9.999925e-01 7.452330e-06
#> [34,]  2.476277e-96 9.805847e-01 1.941534e-02
#> [35,] 4.143575e-104 9.873520e-01 1.264803e-02
#> [36,] 1.259738e-132 5.388373e-01 4.611627e-01
#> [37,] 1.516971e-113 8.572528e-01 1.427472e-01
#> [38,]  5.948781e-79 9.999966e-01 3.360240e-06
#> [39,] 5.064064e-103 9.999307e-01 6.934028e-05
#> [40,]  2.736888e-98 9.939288e-01 6.071194e-03
#> [41,]  2.886594e-56 9.999813e-01 1.867254e-05
#> [42,]  5.138169e-62 9.999875e-01 1.254357e-05
#> [43,]  1.408357e-95 9.524294e-01 4.757059e-02
#> [44,] 4.782959e-106 9.999778e-01 2.221000e-05
#> [45,]  1.234604e-70 9.997790e-01 2.209950e-04
#> [46,]  2.291561e-73 9.793616e-01 2.063839e-02
#> [47,]  2.984550e-65 9.999614e-01 3.863606e-05
#> [48,]  1.512584e-75 9.986823e-01 1.317728e-03
#> [49,]  5.358215e-74 9.995505e-01 4.494678e-04
#> [50,]  1.245663e-39 9.999997e-01 2.512703e-07
#> [51,] 7.729266e-266 6.791296e-09 1.000000e+00
#> [52,] 9.954921e-213 1.509018e-04 9.998491e-01
#> [53,] 1.897537e-218 9.106944e-06 9.999909e-01
#> [54,] 4.128232e-117 5.987163e-04 9.994013e-01
#> [55,] 3.142154e-175 1.044183e-03 9.989558e-01
#> [56,] 9.517474e-167 1.477429e-02 9.852257e-01
#> [57,] 8.912631e-196 1.682074e-04 9.998318e-01
#> [58,] 1.761936e-219 4.989459e-14 1.000000e+00
#> [59,] 4.579693e-159 8.943397e-02 9.105660e-01
#> [60,] 2.716293e-289 2.797982e-09 1.000000e+00
#> [61,] 2.332060e-229 1.500825e-06 9.999985e-01
#> [62,] 7.613591e-240 6.205542e-07 9.999994e-01
#> [63,] 7.145350e-201 4.497101e-03 9.955029e-01
#> [64,] 5.793042e-136 1.159599e-01 8.840401e-01
#> [65,] 6.507162e-199 2.002077e-05 9.999800e-01
#> [66,] 6.845274e-200 1.986525e-04 9.998013e-01
#> [67,] 8.914835e-212 5.129363e-07 9.999995e-01
#> [68,] 2.807218e-125 8.663816e-04 9.991336e-01
#> [69,] 6.459918e-236 1.247735e-07 9.999999e-01
#> [70,] 3.617450e-134 2.216457e-01 7.783543e-01
#> [71,] 3.068336e-239 2.072899e-09 1.000000e+00
#> [72,] 2.825419e-156 6.821447e-04 9.993179e-01
#> [73,] 4.602177e-257 3.672936e-10 1.000000e+00
#> [74,] 9.978143e-156 1.051037e-03 9.989490e-01
#> [75,] 4.631299e-214 2.854840e-06 9.999971e-01
#> 
predTest <- predict(irisMclustDA, X.test)
predTest
#> $classification
#>  [1] setosa     setosa     setosa     setosa     setosa     setosa    
#>  [7] setosa     setosa     setosa     setosa     setosa     setosa    
#> [13] setosa     setosa     setosa     setosa     setosa     setosa    
#> [19] setosa     setosa     setosa     setosa     setosa     setosa    
#> [25] setosa     versicolor versicolor versicolor versicolor versicolor
#> [31] versicolor versicolor versicolor versicolor versicolor versicolor
#> [37] versicolor versicolor versicolor versicolor versicolor virginica 
#> [43] versicolor versicolor versicolor versicolor versicolor versicolor
#> [49] versicolor versicolor virginica  virginica  virginica  virginica 
#> [55] virginica  virginica  virginica  virginica  virginica  virginica 
#> [61] virginica  virginica  virginica  virginica  virginica  virginica 
#> [67] versicolor virginica  virginica  virginica  virginica  virginica 
#> [73] virginica  virginica  virginica 
#> Levels: setosa versicolor virginica
#> 
#> $z
#>              setosa   versicolor    virginica
#>  [1,]  1.000000e+00 2.523334e-18 6.021502e-34
#>  [2,]  1.000000e+00 2.881796e-18 3.538194e-31
#>  [3,]  1.000000e+00 1.178420e-24 4.486508e-40
#>  [4,]  1.000000e+00 2.553955e-22 1.305100e-36
#>  [5,]  1.000000e+00 4.143661e-20 1.376911e-33
#>  [6,]  1.000000e+00 5.893258e-21 2.158059e-33
#>  [7,]  1.000000e+00 1.073792e-16 6.237856e-31
#>  [8,]  1.000000e+00 1.328461e-33 9.307663e-52
#>  [9,]  1.000000e+00 1.288713e-22 7.980067e-39
#> [10,]  1.000000e+00 4.729930e-26 3.717387e-41
#> [11,]  1.000000e+00 5.586897e-22 9.882706e-38
#> [12,]  1.000000e+00 1.162587e-12 1.425626e-27
#> [13,]  1.000000e+00 9.587273e-17 4.827376e-31
#> [14,]  1.000000e+00 6.510286e-25 6.158055e-40
#> [15,]  1.000000e+00 1.467071e-18 7.029886e-31
#> [16,]  1.000000e+00 4.325261e-20 1.852521e-37
#> [17,]  1.000000e+00 6.714034e-37 2.755097e-53
#> [18,]  1.000000e+00 5.423357e-21 2.226046e-38
#> [19,]  1.000000e+00 3.761178e-27 1.810110e-40
#> [20,]  1.000000e+00 7.274004e-23 1.269579e-37
#> [21,]  9.499006e-01 5.009940e-02 9.802674e-19
#> [22,]  1.000000e+00 3.394655e-12 4.209226e-27
#> [23,]  1.000000e+00 7.027520e-15 1.734097e-30
#> [24,]  1.000000e+00 1.376843e-19 3.713174e-33
#> [25,]  1.000000e+00 2.351319e-21 9.432549e-37
#> [26,]  2.680809e-97 9.996813e-01 3.187179e-04
#> [27,]  1.795491e-74 9.992765e-01 7.234846e-04
#> [28,]  7.159623e-82 9.887687e-01 1.123129e-02
#> [29,]  9.103698e-41 9.999938e-01 6.174591e-06
#> [30,]  3.944584e-75 9.917816e-01 8.218400e-03
#> [31,]  1.988148e-89 9.989821e-01 1.017940e-03
#> [32,]  1.833995e-95 9.881890e-01 1.181101e-02
#> [33,]  7.268336e-89 9.999991e-01 9.010643e-07
#> [34,]  5.321069e-53 9.995583e-01 4.417260e-04
#> [35,]  2.461511e-57 9.999653e-01 3.471856e-05
#> [36,]  1.099304e-70 9.999884e-01 1.164699e-05
#> [37,]  3.664552e-80 9.426941e-01 5.730585e-02
#> [38,]  3.819025e-89 9.999965e-01 3.536921e-06
#> [39,] 1.479980e-132 8.652486e-01 1.347514e-01
#> [40,]  2.936411e-42 9.999999e-01 1.179693e-07
#> [41,]  2.136796e-48 9.999949e-01 5.094391e-06
#> [42,] 8.404582e-126 2.987918e-01 7.012082e-01
#> [43,] 3.826070e-104 9.990075e-01 9.925433e-04
#> [44,]  8.407400e-86 9.996890e-01 3.109789e-04
#> [45,]  1.402383e-72 9.993712e-01 6.288146e-04
#> [46,]  4.877834e-92 9.968721e-01 3.127878e-03
#> [47,]  5.372098e-42 9.999979e-01 2.059269e-06
#> [48,]  1.510319e-65 9.996610e-01 3.389848e-04
#> [49,]  1.558161e-77 9.999034e-01 9.656858e-05
#> [50,]  1.464651e-72 9.997568e-01 2.432354e-04
#> [51,] 2.825419e-156 6.821447e-04 9.993179e-01
#> [52,] 2.549967e-162 2.311842e-02 9.768816e-01
#> [53,] 1.400504e-245 5.089164e-06 9.999949e-01
#> [54,] 2.148191e-195 6.863124e-05 9.999314e-01
#> [55,] 1.275329e-274 5.468241e-07 9.999995e-01
#> [56,] 4.818859e-165 3.321108e-03 9.966789e-01
#> [57,] 9.209048e-166 2.096834e-06 9.999979e-01
#> [58,] 2.036128e-212 1.062788e-07 9.999999e-01
#> [59,] 1.025724e-258 6.938810e-03 9.930612e-01
#> [60,] 2.563761e-118 3.305876e-01 6.694124e-01
#> [61,] 2.439902e-160 1.304494e-05 9.999870e-01
#> [62,] 3.870421e-140 7.557603e-02 9.244240e-01
#> [63,] 5.110005e-180 9.441091e-03 9.905589e-01
#> [64,] 1.434938e-137 2.542019e-01 7.457981e-01
#> [65,] 4.969278e-153 5.967287e-03 9.940327e-01
#> [66,] 7.338521e-221 1.461263e-02 9.853874e-01
#> [67,] 6.517202e-117 6.762001e-01 3.237999e-01
#> [68,] 5.549239e-251 1.777694e-07 9.999998e-01
#> [69,] 9.489170e-158 1.638715e-01 8.361285e-01
#> [70,] 7.915383e-192 4.097626e-04 9.995902e-01
#> [71,] 6.964811e-207 1.427441e-08 1.000000e+00
#> [72,] 1.546248e-236 3.756445e-06 9.999962e-01
#> [73,] 1.421391e-210 6.961792e-09 1.000000e+00
#> [74,] 1.791637e-171 2.458719e-03 9.975413e-01
#> [75,] 3.574474e-143 1.192001e-01 8.807999e-01
#> 
# }
```
