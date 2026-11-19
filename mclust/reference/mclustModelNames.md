# MCLUST Model Names

Description of model names used in the *MCLUST* package.

## Usage

``` r
mclustModelNames(model)
```

## Arguments

- model:

  A string specifying the model.

## Details

The following models are available in package mclust:\

**univariate mixture**\

- `"E"`:

  equal variance (one-dimensional)

- `"V"`:

  variable/unqual variance (one-dimensional)

**multivariate mixture**\

- `"EII"`:

  spherical, equal volume

- `"VII"`:

  spherical, unequal volume

- `"EEI"`:

  diagonal, equal volume and shape

- `"VEI"`:

  diagonal, varying volume, equal shape

- `"EVI"`:

  diagonal, equal volume, varying shape

- `"VVI"`:

  diagonal, varying volume and shape

- `"EEE"`:

  ellipsoidal, equal volume, shape, and orientation

- `"VEE"`:

  ellipsoidal, equal shape and orientation (\*)

- `"EVE"`:

  ellipsoidal, equal volume and orientation (\*)

- `"VVE"`:

  ellipsoidal, equal orientation (\*)

- `"EEV"`:

  ellipsoidal, equal volume and equal shape

- `"VEV"`:

  ellipsoidal, equal shape

- `"EVV"`:

  ellipsoidal, equal volume (\*)

- `"VVV"`:

  ellipsoidal, varying volume, shape, and orientation

**single component**\

- `"X"`:

  univariate normal

- `"XII"`:

  spherical multivariate normal

- `"XXI"`:

  diagonal multivariate normal

- `"XXX"`:

  ellipsoidal multivariate normal

(\*) new models in mclust version \>= 5.0.0.

## Value

Returns a list with the following components:

- model:

  a character string indicating the model (as in input).

- type:

  the description of the indicated model (see Details section).

## See also

[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.md),
[`mclustBIC`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)

## Examples

``` r
mclustModelNames("E")
#> $model
#> [1] "E"
#> 
#> $type
#> [1] "univariate, equal variance"
#> 
mclustModelNames("EEE")
#> $model
#> [1] "EEE"
#> 
#> $type
#> [1] "ellipsoidal, equal volume, shape and orientation"
#> 
mclustModelNames("VVV")
#> $model
#> [1] "VVV"
#> 
#> $type
#> [1] "ellipsoidal, varying volume, shape, and orientation"
#> 
mclustModelNames("XXI")
#> $model
#> [1] "XXI"
#> 
#> $type
#> [1] "diagonal multivariate normal"
#> 
```
