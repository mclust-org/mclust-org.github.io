# Plot Combined Clusterings Results

Plot combined clusterings results: classifications corresponding to
`Mclust`/BIC and to the hierarchically combined classes, "entropy plots"
to help to select a number of classes, and the tree structure obtained
from combining mixture components.

## Usage

``` r
# S3 method for class 'clustCombi'
plot(x, what = c("classification", "entropy", "tree"), ...)
```

## Arguments

- x:

  Object returned by
  [`clustCombi`](https://mclust-org.github.io/mclust/reference/clustCombi.md)
  function.

- what:

  Type of plot.

- ...:

  Other arguments to be passed to other functions:
  [`combiPlot`](https://mclust-org.github.io/mclust/reference/combiPlot.md),
  [`entPlot`](https://mclust-org.github.io/mclust/reference/entPlot.md),
  [`combiTree`](https://mclust-org.github.io/mclust/reference/combiTree.md).
  Please see the corresponding documentations.

## Value

Classifications are plotted with
[`combiPlot`](https://mclust-org.github.io/mclust/reference/combiPlot.md),
which relies on the `Mclust` plot functions. Entropy plots are plotted
with
[`entPlot`](https://mclust-org.github.io/mclust/reference/entPlot.md)
and may help to select a number of classes: please see the article cited
in the references. Tree plots are produced by
[`combiTree`](https://mclust-org.github.io/mclust/reference/combiTree.md)
and graph the tree structure implied by the clusters combining process.

## References

J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010).
Combining mixture components for clustering. *Journal of Computational
and Graphical Statistics, 19(2):332-353.*

## Author

J.-P. Baudry, A. E. Raftery, L. Scrucca

## See also

[`combiPlot`](https://mclust-org.github.io/mclust/reference/combiPlot.md),
[`entPlot`](https://mclust-org.github.io/mclust/reference/entPlot.md),
[`combiTree`](https://mclust-org.github.io/mclust/reference/combiTree.md),
[`clustCombi`](https://mclust-org.github.io/mclust/reference/clustCombi.md).

## Examples

``` r
# \donttest{
data(Baudry_etal_2010_JCGS_examples)

## 1D Example 
output <- clustCombi(data = Test1D, G=1:15)

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output) 














## 2D Example 
output <- clustCombi(data = ex4.1) 

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output) 










## 3D Example 
output <- clustCombi(data = ex4.4.2)

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes (please see the article cited 
# in the references)
plot(output)






# }
```
