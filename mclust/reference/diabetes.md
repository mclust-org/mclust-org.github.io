# Diabetes Data (flawed)

The data set contains three measurements made on 145 non-obese adult
patients classified into three groups.

## Usage

``` r
data(diabetes)
```

## Format

A data frame with the following variables:

- class:

  The type of diabete: `Normal`, `Overt`, and `Chemical`.

- glucose:

  Area under plasma glucose curve after a three hour oral glucose
  tolerance test (OGTT).

- insulin:

  Area under plasma insulin curve after a three hour oral glucose
  tolerance test (OGTT).

- sspg:

  Steady state plasma glucose.

## Details

This dataset is *flawed* (compare with the reference) and it is provided
here only for backward compatibility. A 5-variable version of the Reaven
and Miller data is available in package rrcov. The *glucose* and *sspg*
columns in this dataset are identical to the *fpg* and *insulin*
columns, respectively in the rrcov version. The *insulin* column in this
dataset differs from the *glucose* column in the rrcov version in one
entry: observation 104 has the value 45 in the *insulin* column in this
data, and 455 in the corresponding *glucose* column of the rrcov
version.

## Source

Reaven, G. M. and Miller, R. G. (1979). An attempt to define the nature
of chemical diabetes using a multidimensional analysis. *Diabetologia*
16:17-24.
