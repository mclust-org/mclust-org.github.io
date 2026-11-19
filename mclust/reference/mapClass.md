# Correspondence between classifications

Best correspondence between classes given two vectors viewed as
alternative classifications of the same object.

## Usage

``` r
mapClass(a, b)
```

## Arguments

- a:

  A numeric or character vector of class labels.

- b:

  A numeric or character vector of class labels. Must have the same
  length as `a`.

## Value

A list with two named elements, `aTOb` and `bTOa` which are themselves
lists. The `aTOb` list has a component corresponding to each unique
element of `a`, which gives the element or elements of `b` that result
in the closest class correspondence.

The `bTOa` list has a component corresponding to each unique element of
`b`, which gives the element or elements of `a` that result in the
closest class correspondence.

## See also

[`classError`](https://mclust-org.github.io/mclust/reference/classError.md),
[`table`](https://rdrr.io/r/base/table.html)

## Examples

``` r
a <- rep(1:3, 3)
a
#> [1] 1 2 3 1 2 3 1 2 3
b <- rep(c("A", "B", "C"), 3)
b
#> [1] "A" "B" "C" "A" "B" "C" "A" "B" "C"
mapClass(a, b)
#> $aTOb
#> $aTOb$`1`
#> [1] "A"
#> 
#> $aTOb$`2`
#> [1] "B"
#> 
#> $aTOb$`3`
#> [1] "C"
#> 
#> 
#> $bTOa
#> $bTOa$A
#> [1] 1
#> 
#> $bTOa$B
#> [1] 2
#> 
#> $bTOa$C
#> [1] 3
#> 
#> 
a <- sample(1:3, 9, replace = TRUE)
a
#> [1] 2 3 2 3 2 1 3 3 1
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
#> [1] "C" "A" "A" "C" "B" "C" "B" "A" "B"
mapClass(a, b)
#> $aTOb
#> $aTOb$`1`
#> [1] "B"
#> 
#> $aTOb$`2`
#> [1] "A"
#> 
#> $aTOb$`3`
#> [1] "A"
#> 
#> 
#> $bTOa
#> $bTOa$A
#> [1] 3
#> 
#> $bTOa$B
#> [1] 1
#> 
#> $bTOa$C
#> [1] 1
#> 
#> 
```
