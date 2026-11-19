# Majority vote

A function to compute the majority vote (some would say plurality) label
in a vector of labels, breaking ties at random.

## Usage

``` r
majorityVote(x)
```

## Arguments

- x:

  A vector of values, either numerical or not.

## Value

A list with the following components:

- table:

  A table of votes for each unique value of `x`.

- ind:

  An integer specifying which unique value of `x` corresponds to the
  majority vote.

- majority:

  A string specifying the majority vote label.

## Author

L. Scrucca

## Examples

``` r
x <- c("A", "C", "A", "B", "C", "B", "A")
majorityVote(x)
#> $table
#> x
#> A B C 
#> 3 2 2 
#> 
#> $ind
#> [1] 1
#> 
#> $majority
#> [1] "A"
#> 
```
