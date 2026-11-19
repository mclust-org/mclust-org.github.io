# Log-likelihood from a table of BIC values for parameterized Gaussian mixture models

Compute the maximal log-likelihood from a table of BIC values contained
in a `'mclustBIC'` object as returned by function
[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).

## Usage

``` r
mclustLoglik(object, ...)
```

## Arguments

- object:

  An object of class `'mclustBIC'` containing the BIC values as returned
  by a call to
  [`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).

- ...:

  Catches unused arguments in an indirect or list call via `do.call`.

## Value

An object of class `'mclustLoglik'` containing the maximal
log-likelihood values for the Gaussian mixture models provided as input.

## See also

[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).

## Examples

``` r
# \donttest{
BIC <- mclustBIC(iris[,1:4])
mclustLoglik(BIC)
#> Log-likelihood: 
#>         EII       VII       EEI       VEI       EVI       VVI       EEE
#> 1 -889.5161 -889.5161 -741.0175 -741.0175 -741.0175 -741.0175 -379.9146
#> 2 -536.6527 -478.5591 -488.9148 -443.0667 -463.5690 -386.1853 -296.4476
#> 3 -401.8027 -384.3168 -361.4295 -339.4719 -338.7895 -307.1808 -256.3547
#> 4 -396.7007 -348.6801 -356.0795 -309.0882 -338.6024 -287.8238 -250.3587
#> 5 -328.6891 -298.6500 -300.8103 -264.0030 -283.1952 -245.4911 -217.2257
#> 6 -282.6973 -265.2044 -264.2199 -242.8828 -266.7784 -220.8632 -207.2198
#> 7 -278.2550 -246.5526 -261.3891 -230.1349 -266.4632 -227.9953 -206.0134
#> 8 -262.3275 -232.9904 -237.9780 -214.4661 -209.8631 -204.2209 -196.8715
#> 9 -234.5210 -217.2319 -227.8751 -210.7094 -203.5193 -177.4891 -191.2567
#>         VEE       EVE       VVE       EEV       VEV       EVV       VVV
#> 1 -379.9146 -379.9146 -379.9146 -379.9146 -379.9146 -379.9146 -379.9146
#> 2 -278.0572 -273.4962 -244.9697 -259.6669 -215.7260 -259.0164 -214.3547
#> 3 -237.5609 -258.1150 -238.0428 -232.1991 -186.0740 -222.7946 -180.1858
#> 4 -222.2484 -257.5697 -217.1359 -232.1843 -175.7393 -222.3485 -167.4862
#> 5        NA -246.6153 -190.8369 -180.8395 -161.8148        NA -152.9095
#> 6 -194.5452 -195.6877 -184.5974 -159.3600 -154.3623        NA -154.4236
#> 7        NA -194.4254 -174.9438 -154.5511 -136.9276 -159.3927 -142.9108
#> 8 -187.1141 -174.8416 -162.0602 -139.7393 -110.9183 -135.2804 -117.1858
#> 9        NA -169.4178 -151.6925 -123.9250  -98.8269 -125.5495 -106.1340
# }
```
