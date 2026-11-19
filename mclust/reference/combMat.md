# Combining Matrix

Create a combining matrix

## Usage

``` r
combMat(K, l1, l2)
```

## Arguments

- K:

  The original number of classes: the matrix will define a combining
  from K to (K-1) classes.

- l1:

  Label of one of the two classes to be combined.

- l2:

  Label of the other class to be combined.

## Value

If `z` is a vector (length *K*) whose *k*th entry is the probability
that an observation belongs to the *k*th class in a *K*-classes
classification, then `combiM %*% z` is the vector (length *K-1*) whose
*k*th entry is the probability that the observation belongs to the *k*th
class in the *K-1*-classes classification obtained by merging classes
`l1` and `l2` in the initial classification.

## Author

J.-P. Baudry, A. E. Raftery, L. Scrucca

## See also

[`clustCombi`](https://mclust-org.github.io/mclust/reference/clustCombi.md),
[`combiPlot`](https://mclust-org.github.io/mclust/reference/combiPlot.md)
