# UCI Wisconsin Diagnostic Breast Cancer Data

The data set provides data for 569 patients on 30 features of the cell
nuclei obtained from a digitized image of a fine needle aspirate (FNA)
of a breast mass. For each patient the cancer was diagnosed as malignant
or benign.

## Usage

``` r
data(wdbc)
```

## Format

A data frame with 569 observations on the following variables:

- `ID`:

  ID number

- `Diagnosis`:

  cancer diagnosis: `M` = malignant, `B` = benign

- `Radius_mean`:

  a numeric vector

- `Texture_mean`:

  a numeric vector

- `Perimeter_mean`:

  a numeric vector

- `Area_mean`:

  a numeric vector

- `Smoothness_mean`:

  a numeric vector

- `Compactness_mean`:

  a numeric vector

- `Concavity_mean`:

  a numeric vector

- `Nconcave_mean`:

  a numeric vector

- `Symmetry_mean`:

  a numeric vector

- `Fractaldim_mean`:

  a numeric vector

- `Radius_se`:

  a numeric vector

- `Texture_se`:

  a numeric vector

- `Perimeter_se`:

  a numeric vector

- `Area_se`:

  a numeric vector

- `Smoothness_se`:

  a numeric vector

- `Compactness_se`:

  a numeric vector

- `Concavity_se`:

  a numeric vector

- `Nconcave_se`:

  a numeric vector

- `Symmetry_se`:

  a numeric vector

- `Fractaldim_se`:

  a numeric vector

- `Radius_extreme`:

  a numeric vector

- `Texture_extreme`:

  a numeric vector

- `Perimeter_extreme`:

  a numeric vector

- `Area_extreme`:

  a numeric vector

- `Smoothness_extreme`:

  a numeric vector

- `Compactness_extreme`:

  a numeric vector

- `Concavity_extreme`:

  a numeric vector

- `Nconcave_extreme`:

  a numeric vector

- `Symmetry_extreme`:

  a numeric vector

- `Fractaldim_extreme`:

  a numeric vector

## Details

The recorded features are:

- `Radius` as mean of distances from center to points on the perimeter

- `Texture` as standard deviation of gray-scale values

- `Perimeter` as cell nucleus perimeter

- `Area` as cell nucleus area

- `Smoothness` as local variation in radius lengths

- `Compactness` as cell nucleus compactness, perimeter^2 / area - 1

- `Concavity` as severity of concave portions of the contour

- `Nconcave` as number of concave portions of the contour

- `Symmetry` as cell nucleus shape

- `Fractaldim` as fractal dimension, "coastline approximation" - 1

For each feature the recorded values are computed from each image as
`<feature_name>_mean`, `<feature_name>_se`, and
`<feature_name>_extreme`, for the mean, the standard error, and the mean
of the three largest values.

## Source

The Breast Cancer Wisconsin (Diagnostic) Data Set (`wdbc.data`,
`wdbc.names`) from the UCI Machine Learning Repository
<https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)>.
Please note the UCI conditions of use.

## References

Mangasarian, O. L., Street, W. N., and Wolberg, W. H. (1995) Breast
cancer diagnosis and prognosis via linear programming. *Operations
Research*, 43(4), pp. 570-577.
