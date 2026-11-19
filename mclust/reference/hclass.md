# Classifications from Hierarchical Agglomeration

Determines the classifications corresponding to different numbers of
groups given merge pairs from hierarchical agglomeration.

## Usage

``` r
hclass(hcPairs, G)
```

## Arguments

- hcPairs:

  A numeric two-column matrix in which the *i*th row gives the minimum
  index for observations in each of the two clusters merged at the *i*th
  stage of agglomerative hierarchical clustering.

- G:

  An integer or vector of integers giving the number of clusters for
  which the corresponding classfications are wanted.

## Value

A matrix with `length(G)` columns, each column corresponding to a
classification. Columns are indexed by the character representation of
the integers in `G`.

## See also

[`hc`](https://mclust-org.github.io/mclust/reference/hc.md),
[`hcE`](https://mclust-org.github.io/mclust/reference/hcE.md)

## Examples

``` r
hcTree <- hc(modelName="VVV", data = iris[,-5])
cl <- hclass(hcTree,c(2,3))

# \donttest{
par(pty = "s", mfrow = c(1,1))
clPairs(iris[,-5],cl=cl[,"2"])

clPairs(iris[,-5],cl=cl[,"3"])

# }
```
