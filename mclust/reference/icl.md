# ICL for an estimated Gaussian Mixture Model

Computes the ICL (Integrated Complete-data Likelihood) for criterion for
a Gaussian Mixture Model fitted by
[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md).

## Usage

``` r
icl(object, ...)
```

## Arguments

- object:

  An object of class `'Mclust'` resulting from a call to
  [`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md).

- ...:

  Further arguments passed to or from other methods.

## Value

The ICL for the given input MCLUST model.

## References

Biernacki, C., Celeux, G., Govaert, G. (2000). Assessing a mixture model
for clustering with the integrated completed likelihood. *IEEE Trans.
Pattern Analysis and Machine Intelligence*, 22 (7), 719-725.

## See also

[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md),
[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md),
[`mclustICL`](https://mclust-org.github.io/mclust/reference/mclustICL.md),
[`bic`](https://mclust-org.github.io/mclust/reference/bic.md).

## Examples

``` r
mod <- Mclust(iris[,1:4])
icl(mod)
#> [1] -561.7289
```
