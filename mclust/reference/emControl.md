# Set control values for use with the EM algorithm

Supplies a list of values including tolerances for singularity and
convergence assessment, for use functions involving EM within *MCLUST*.

## Usage

``` r
emControl(eps, tol, itmax, equalPro)
```

## Arguments

- eps:

  A scalar tolerance associated with deciding when to terminate
  computations due to computational singularity in covariances. Smaller
  values of `eps` allow computations to proceed nearer to singularity.
  The default is the relative machine precision `.Machine$double.eps`,
  which is approximately \\(2e-16\\) on IEEE-compliant machines.

- tol:

  A vector of length two giving relative convergence tolerances for the
  log-likelihood and for parameter convergence in the inner loop for
  models with iterative M-step ("VEI", "VEE", "EVE", "VVE", "VEV"),
  respectively. The default is `c(1.e-5, sqrt(.Machine$double.eps))`. If
  only one number is supplied, it is used as the tolerance for the outer
  iterations and the tolerance for the inner iterations is as in the
  default.

- itmax:

  A vector of length two giving integer limits on the number of EM
  iterations and on the number of iterations in the inner loop for
  models with iterative M-step ("VEI", "VEE", "EVE", "VVE", "VEV"),
  respectively. The default is
  `c(.Machine$integer.max, .Machine$integer.max)` allowing termination
  to be completely governed by `tol`. If only one number is supplied, it
  is used as the iteration limit for the outer iteration only.

- equalPro:

  Logical variable indicating whether or not the mixing proportions are
  equal in the model. Default: `equalPro = FALSE`.

## Value

A named list in which the names are the names of the arguments and the
values are the values supplied to the arguments.

## Details

`emControl` is provided for assigning values and defaults for EM within
*MCLUST*.

## See also

[`em`](https://mclust-org.github.io/mclust/reference/em.md),
[`estep`](https://mclust-org.github.io/mclust/reference/estep.md),
[`me`](https://mclust-org.github.io/mclust/reference/me.md),
[`mstep`](https://mclust-org.github.io/mclust/reference/mstep.md),
[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)

## Examples

``` r
irisBIC <- mclustBIC(iris[,-5], control = emControl(tol = 1.e-6))
summary(irisBIC, iris[,-5])
#> Best BIC values:
#>              VEV,2        VEV,3      VVV,2
#> BIC      -561.7285 -562.5508708 -574.01783
#> BIC diff    0.0000   -0.8224087  -12.28937
#> 
#> Classification table for model (VEV,2): 
#> 
#>   1   2 
#>  50 100 
```
