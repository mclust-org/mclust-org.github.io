# UCI Thyroid Gland Data

Data on five laboratory tests administered to a sample of 215 patients.
The tests are used to predict whether a patient's thyroid can be
classified as euthyroidism (normal thyroid gland function),
hypothyroidism (underactive thyroid not producing enough thyroid
hormone) or hyperthyroidism (overactive thyroid producing and secreting
excessive amounts of the free thyroid hormones T3 and/or thyroxine T4).
Diagnosis of thyroid operation was based on a complete medical record,
including anamnesis, scan, etc.

## Usage

``` r
data(thyroid)
```

## Format

A data frame with the following variables:

- Diagnosis:

  Diagnosis of thyroid operation: `Hypo`, `Normal`, and `Hyper`.

- RT3U:

  T3-resin uptake test (percentage).

- T4:

  Total Serum thyroxin as measured by the isotopic displacement method.

- T3:

  Total serum triiodothyronine as measured by radioimmuno assay.

- TSH:

  Basal thyroid-stimulating hormone (TSH) as measured by radioimmuno
  assay.

- DTSH:

  Maximal absolute difference of TSH value after injection of 200 micro
  grams of thyrotropin-releasing hormone as compared to the basal value.

## Source

One of several databases in the Thyroid Disease Data Set
(`new-thyroid.data`, `new-thyroid.names`) of the UCI Machine Learning
Repository <https://archive.ics.uci.edu/ml/datasets/thyroid+disease>.
Please note the UCI conditions of use.

## References

Coomans, D., Broeckaert, M. Jonckheer M. and Massart D.L. (1983)
Comparison of Multivariate Discriminant Techniques for Clinical Data -
Application to the Thyroid Functional State, *Meth. Inform. Med.* 22,
pp. 93-101.

Coomans, D. and I. Broeckaert (1986) *Potential Pattern Recognition in
Cemical and Medical Decision Making*, Research Studies Press,
Letchworth, England.
