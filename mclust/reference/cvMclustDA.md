# MclustDA cross-validation

V-fold cross-validation for classification models based on Gaussian
finite mixture modelling.

## Usage

``` r
cvMclustDA(object, nfold = 10, 
           prop = object$prop,
           verbose = interactive(), 
           ...)
```

## Arguments

- object:

  An object of class `'MclustDA'` resulting from a call to
  [`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).

- nfold:

  An integer specifying the number of folds (by defaul 10-fold CV is
  used).

- prop:

  A vector of class prior probabilities, which if not provided default
  to the class proportions in the training data.

- verbose:

  A logical controlling if a text progress bar is displayed during the
  cross-validation procedure. By default is `TRUE` if the session is
  interactive, and `FALSE` otherwise.

- ... :

  Further arguments passed to or from other methods.

## Details

The function implements V-fold cross-validation for classification
models fitted by
[`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md).
Classification error and Brier score are the metrics returned, but other
metrics can be computed using the output returned by this function (see
Examples section below).

## Value

The function returns a list with the following components:

- classification:

  a factor of cross-validated class labels.

- z:

  a matrix containing the cross-validated probabilites for class
  assignment.

- ce:

  the cross-validation classification error.

- se.ce:

  the standard error of the cross-validated classification error.

- brier:

  the cross-validation Brier score.

- se.brier:

  the standard error of the cross-validated Brier score.

## Author

Luca Scrucca

## See also

[`MclustDA`](https://mclust-org.github.io/mclust/reference/MclustDA.md),
[`predict.MclustDA`](https://mclust-org.github.io/mclust/reference/predict.MclustDA.md),
[`classError`](https://mclust-org.github.io/mclust/reference/classError.md),
[`BrierScore`](https://mclust-org.github.io/mclust/reference/BrierScore.md)

## Examples

``` r
# \donttest{
# Iris data
Class <- iris$Species
X <- iris[,1:4]

## EDDA model with common covariance (essentially equivalent to linear discriminant analysis)
irisEDDA <- MclustDA(X, Class, modelType = "EDDA", modelNames = "EEE")
cv <- cvMclustDA(irisEDDA)                         # 10-fold CV (default)
str(cv)
#> List of 6
#>  $ classification: Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ z             : num [1:150, 1:3] 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "setosa" "versicolor" "virginica"
#>  $ ce            : num 0.02
#>  $ se.ce         : num 0.0102
#>  $ brier         : num 0.0187
#>  $ se.brier      : num 0.00703
cv <- cvMclustDA(irisEDDA, nfold = length(Class))  # LOO-CV
str(cv)
#> List of 6
#>  $ classification: Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ z             : num [1:150, 1:3] 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "setosa" "versicolor" "virginica"
#>  $ ce            : num 0.02
#>  $ se.ce         : num 0.0115
#>  $ brier         : num 0.017
#>  $ se.brier      : num 0.00824

## MclustDA model selected by BIC
irisMclustDA <- MclustDA(X, Class)
cv <- cvMclustDA(irisMclustDA)                     # 10-fold CV (default)
str(cv)
#> List of 6
#>  $ classification: Factor w/ 3 levels "setosa","versicolor",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ z             : num [1:150, 1:3] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "setosa" "versicolor" "virginica"
#>  $ ce            : num 0.02
#>  $ se.ce         : num 0.0102
#>  $ brier         : num 0.0174
#>  $ se.brier      : num 0.00712

# Banknote data
data("banknote")
Class <- banknote$Status
X <- banknote[,2:7]

## EDDA model selected by BIC
banknoteEDDA <- MclustDA(X, Class, modelType = "EDDA")
cv <- cvMclustDA(banknoteEDDA)                     # 10-fold CV (default)
str(cv)
#> List of 6
#>  $ classification: Factor w/ 2 levels "counterfeit",..: 2 2 2 2 2 2 2 2 2 2 ...
#>  $ z             : num [1:200, 1:2] 2.45e-06 5.47e-25 4.80e-23 1.20e-24 7.38e-21 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "counterfeit" "genuine"
#>  $ ce            : num 0.005
#>  $ se.ce         : num 0.005
#>  $ brier         : num 0.00499
#>  $ se.brier      : num 0.00494

(ConfusionMatrix <- table(Pred = cv$classification, Class))
#>              Class
#> Pred          counterfeit genuine
#>   counterfeit         100       1
#>   genuine               0      99
TP <- ConfusionMatrix[1,1]
FP <- ConfusionMatrix[1,2]
FN <- ConfusionMatrix[2,1]
TN <- ConfusionMatrix[2,2]
(Sensitivity <- TP/(TP+FN))
#> [1] 1
(Specificity <- TN/(FP+TN))
#> [1] 0.99
# }
```
