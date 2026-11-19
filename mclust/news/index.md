# Changelog

## mclust 6.1.4

CRAN release: 2026-11-18

- Fixed issues with structure() formatting as reported by CRAN.

## mclust 6.1.3

CRAN release: 2026-07-05

- Bug fix on
  [`mclustBootstrapLRT()`](https://mclust-org.github.io/mclust/reference/mclustBootstrapLRT.md)
  when model under null hypothesis is XXI or XXI with one component.
- Added a check for dropping constant columns in input data matrix.
- Added [`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
  `Mclust` and `densityMclust` methods.
- Small fixes to mclust.f fortran code to avoid warnings.
- Add count() function to efficiently tabulate frequencies given bins.
- Fix bugs in colnames for 1D case.
- Fix (marginal) issues in some `em*()` functions that end with the
  M-step instead of the E-step.
- Bug fix on
  [`summary.crimcoords()`](https://mclust-org.github.io/mclust/reference/crimcoords.md)
  when `numdir` is provided in
  [`crimcoords()`](https://mclust-org.github.io/mclust/reference/crimcoords.md)
  call.
- Set varnames for input data if missing from the name of the input
  vector or data matrix in all main functions, namely
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md),
  [`densityMclust()`](https://mclust-org.github.io/mclust/reference/densityMclust.md),
  [`MclustDA()`](https://mclust-org.github.io/mclust/reference/MclustDA.md),
  and
  [`MclustSSC()`](https://mclust-org.github.io/mclust/reference/MclustSSC.md).

## mclust 6.1.1

CRAN release: 2024-04-29

- Corrected computation of df for MclustDA and EDDA classification
  models.
- Bug fix on
  [`sim()`](https://mclust-org.github.io/mclust/reference/sim.md) when
  EVV, EVE, and VVE models for G = 1.
- Bug fix on
  [`summary.MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/summary.MclustBootstrap.md)
  when computing confidence intervals for G = 1.

## mclust 6.1

CRAN release: 2024-02-23

- Added
  [`logsumexp()`](https://mclust-org.github.io/mclust/reference/logsumexp.md)
  and
  [`softmax()`](https://mclust-org.github.io/mclust/reference/softmax.md)
  functions as a wrapper to efficiently implementations written in
  Fortran code.
- Substituted R code in several parts of the package with the above
  Fortran-based functions to compute densities and posterior
  probabilities from Gaussian mixtures.
- Changes on
  [`me.weighted()`](https://mclust-org.github.io/mclust/reference/me.weighted.md)
  to use convergence criterion as in other mclust functions and improved
  efficiency by using the above mentioned Fortran-based functions. This
  also brings computational improvements in the weighted likelihood
  bootstrap implemented in `MclustBootstrap(..., type = "wlbs")`.
- Bug fix on
  [`MclustDA()`](https://mclust-org.github.io/mclust/reference/MclustDA.md)
  when number of obs is less than number of vars.
- Bug fix on `MclustBootstrap(object, ..., type = "pb")` when `object`
  is of class `densityMclust`.

## mclust 6.0.1

CRAN release: 2023-11-15

- Changed initialization in
  [`MclustSSC()`](https://mclust-org.github.io/mclust/reference/MclustSSC.md)
  for components of unlabeled data via k-means.
- Corrected output of `summary.MclustSCC()` for components of unlabeled
  data.
- Updated citation info with reference to book published by Chapman &
  Hall/CRC

## mclust 6.0.0

CRAN release: 2022-10-31

- Major release of mclust accompanying the upcoming book by Chapman &
  Hall/CRC.

## mclust 5.4.11

- Added
  [`summary.crimcoords()`](https://mclust-org.github.io/mclust/reference/crimcoords.md)
  method and removed argument `plot` from
  [`crimcoords()`](https://mclust-org.github.io/mclust/reference/crimcoords.md)
  function call.

## mclust 5.4.10

CRAN release: 2022-05-20

- Updated banner on startup.
- Updated info on man page for datasets `diabetes`, `wdbc`, and
  `thyroid`.\
- Std. error for cross-validation in
  [`cvMclustDA()`](https://mclust-org.github.io/mclust/reference/cvMclustDA.md)
  uses formula for the weighted standard deviation with weights given by
  folds size.
- Fix .Rd files.

## mclust 5.4.9

CRAN release: 2021-12-17

- Added
  [`crimcoords()`](https://mclust-org.github.io/mclust/reference/crimcoords.md)
  to compute discriminant coordinates or crimcoords.
- Fixed man page for
  [`cvMclustDA()`](https://mclust-org.github.io/mclust/reference/cvMclustDA.md).

## mclust 5.4.8

CRAN release: 2021-11-05

- [`densityMclust()`](https://mclust-org.github.io/mclust/reference/densityMclust.md)
  by default draw a graph of the density estimate.
- Fixed a bug in computing mixture density if the noise component is
  present.
- Changed default behaviour of
  [`hc()`](https://mclust-org.github.io/mclust/reference/hc.md) when
  called to perform agglomerative hierarchical clustering instead of
  using for EM initialization.
- The default `mclust.options("hcModelName")` now returns only the model
  to be used.
- Changed default `partition` argument of
  [`hc()`](https://mclust-org.github.io/mclust/reference/hc.md) function
  by adding `dupPartion()` to remove data duplicates.
- Added checks to
  [`mclustBootstrapLRT()`](https://mclust-org.github.io/mclust/reference/mclustBootstrapLRT.md)
  to stop if an invalid `modelName` is provided or a one-component
  mixture model is provided.
- Extended the functionality of
  [`cvMclustDA()`](https://mclust-org.github.io/mclust/reference/cvMclustDA.md)
  by including as cross-validated metrics both the classification error
  and the Brier score.
- Updated info on dataset man pages.

## mclust 5.4.7

CRAN release: 2020-11-20

- Updated plot method (dendrogram) for hierarchical clustering — now
  based on classification likelihood.
- Added
  [`MclustSSC()`](https://mclust-org.github.io/mclust/reference/MclustSSC.md)
  function (and related `print`, `summary`, `plot`, and `predict`,
  methods) for semi-supervised classification.
- Exchanged order of models VEE and EVE to account for increasing
  complexity of EVE.
- Added `cex` argument to
  [`clPairs()`](https://mclust-org.github.io/mclust/reference/clPairs.md)
  to control character expansion used in plotting symbols.
- [`em()`](https://mclust-org.github.io/mclust/reference/em.md) and
  [`me()`](https://mclust-org.github.io/mclust/reference/me.md) have now
  `data` as first argument.

## mclust 5.4.6

CRAN release: 2020-04-11

- Fixed issues with source Fortran code with gfortran 10 as reported by
  CRAN.
- Clean code of `hcCriterion()`.
- Replaced `CEX` argument in functions with standard base graph `cex`
  argument.
- Removed `ylim` argument in function; it can be passed via `...`.
- MclustDA models use the default SVD transformation of the data for
  initialization of the EM algorithm.
- Added `icl` criterion to object returned by
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md).
- Fixed number of pages for the RJ reference.
- [`quantileMclust()`](https://mclust-org.github.io/mclust/reference/cdfMclust.md)
  uses bisection line search method for numerically computing quantiles.

## mclust 5.4.5

CRAN release: 2019-07-08

- Fixed warnings in Fortran calls raised by CRAN.

## mclust 5.4.4

CRAN release: 2019-06-27

- Added
  [`classPriorProbs()`](https://mclust-org.github.io/mclust/reference/classPriorProbs.md)
  to estimate prior class probabilities.
- Added
  [`BrierScore()`](https://mclust-org.github.io/mclust/reference/BrierScore.md)
  to compute the Brier score for assessing the accuracy of probabilistic
  predictions.
- Added
  [`randomOrthogonalMatrix()`](https://mclust-org.github.io/mclust/reference/randomOrthogonalMatrix.md)
  to generate random orthogonal basis matrices.
- Partial rewriting of
  [`summary.MclustDA()`](https://mclust-org.github.io/mclust/reference/summary.MclustDA.md)
  internals to provide both the classification error and the Brier score
  for training and/or test data.
- Partial rewriting of
  [`plot.MclustDA()`](https://mclust-org.github.io/mclust/reference/plot.MclustDA.md)
  internals.
- Added
  [`dmvnorm()`](https://mclust-org.github.io/mclust/reference/dmvnorm.md)
  for computing the density of a general multivariate Gaussian
  distribution via efficient Fortran code.
- Added Wisconsin diagnostic breast cancer (WDBC) data.
- Added EuroUnemployment data.
- Fixed mismatches in Fortran calls.
- Bugs fix.

## mclust 5.4.3

CRAN release: 2019-03-14

- Added website site and update DESCRIPTION with URL.
- Fixed a bug when checking for univariate data with a single
  observation in several instances. Using
  [`NCOL()`](https://rdrr.io/r/base/nrow.html) works both for n-values
  vector or nx1 matrix.
- Fixed a bug when `hcPairs` are provided in the `initialization`
  argument of
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)
  (and relatives) and the number of observations exceed the threshold
  for subsetting.
- Fixed bugs on axes for some manual pairs plots.
- Renamed `type = "level"` to `type = "hdr"`, and `level.prob` to
  `prob`, in
  [`surfacePlot()`](https://mclust-org.github.io/mclust/reference/surfacePlot.md)
  for getting HDRs graphs
- Fixed a bug in `type = "hdr"` plot on
  [`surfacePlot()`](https://mclust-org.github.io/mclust/reference/surfacePlot.md).
- Fixed a bug in
  [`as.Mclust()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md).
- Small changes to
  [`summary.MclustDA()`](https://mclust-org.github.io/mclust/reference/summary.MclustDA.md)
  when `modelType = "EDDA"` and in general for a more compact output.

## mclust 5.4.2

CRAN release: 2018-11-17

- Added
  [`mclustBICupdate()`](https://mclust-org.github.io/mclust/reference/mclustBICupdate.md)
  to merge the best values from two BIC results as returned by
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).
- Added
  [`mclustLoglik()`](https://mclust-org.github.io/mclust/reference/mclustLoglik.md)
  to compute the maximal log-likelihood values from BIC results as
  returned by
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).
- Added option `type = "level"` to
  [`plot.densityMclust()`](https://mclust-org.github.io/mclust/reference/plot.densityMclust.md)
  and
  [`surfacePlot()`](https://mclust-org.github.io/mclust/reference/surfacePlot.md)
  to draw highest density regions.
- Added
  [`meXXI()`](https://mclust-org.github.io/mclust/reference/meE.md) and
  [`meXXX()`](https://mclust-org.github.io/mclust/reference/meE.md) to
  exported functions.
- Updated vignette.

## mclust 5.4.1

CRAN release: 2018-06-27

- Added parametric bootstrap option (`type = "pb"`) in
  [`MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md).
- Added the options to get averages of resampling distributions in
  [`summary.MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/summary.MclustBootstrap.md)
  and to plot resampling-based confidence intervals in
  [`plot.MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/plot.MclustBoostrap.md).
- Added function `catwrap()` for wrapping printed lines at
  `getOption("width")` when using
  [`cat()`](https://rdrr.io/r/base/cat.html).
- [`mclust.options()`](https://mclust-org.github.io/mclust/reference/mclust.options.md)
  now modify the variable `.mclust` in the namespace of the package, so
  it should work even inside an mclust-function call.
- Fixed a bug in
  [`covw()`](https://mclust-org.github.io/mclust/reference/covw.md) when
  `normalize = TRUE`.
- Fixed a bug in
  [`estepVEV()`](https://mclust-org.github.io/mclust/reference/estepE.md)
  and
  [`estepVEE()`](https://mclust-org.github.io/mclust/reference/estepE.md)
  when parameters contains `Vinv`.
- Fixed a bug in
  [`plotDensityMclustd()`](https://mclust-org.github.io/mclust/reference/plot.densityMclust.md)
  when drawing marginal axes.
- Fixed a bug in
  [`summary.MclustDA()`](https://mclust-org.github.io/mclust/reference/summary.MclustDA.md)
  when computing classification error in the extreme case of a minor
  class of assignment.
- Fixed a bug in the initialisation of
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)
  when a noise component is present for 1-dimensional data.
- Fixed bugs in some examples documenting
  [`clustCombi()`](https://mclust-org.github.io/mclust/reference/clustCombi.md)
  and related functions.

## mclust 5.4

CRAN release: 2017-11-22

- Model-based hierarchical clustering used to start the EM-algorithm is
  now based on the scaled SVD transformation proposed by Scrucca and
  Raftery (2016). This change is not backward compatible. However,
  previous results can be easily obtained by issuing the command:
  `mclust.options(hcUse = "VARS")` For more details see
  [`help("mclust.options")`](https://mclust-org.github.io/mclust/reference/mclust.options.md).
- Added `subset` parameter in
  [`mclust.options()`](https://mclust-org.github.io/mclust/reference/mclust.options.md)
  to control the maximal sample size to be used in the initial
  model-based hierarchical phase.
- [`predict.densityMclust()`](https://mclust-org.github.io/mclust/reference/predict.densityMclust.md)
  can optionally returns the density on a logarithm scale.
- Removed normalization of mixing proportions for new models in single
  mstep.
- Internal rewrite of code used by
  [`packageStartupMessage()`](https://rdrr.io/r/base/message.html).
- Fixed a small bug in
  [`MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md)
  in the univariate data case.
- Fixed bugs when both the noise and subset are provided for
  initialization.
- Vignette updated to include references, startup message, css style,
  etc.
- Various bug fixes in plotting methods when noise is present.
- Updated references in
  [`citation()`](https://rdrr.io/r/utils/citation.html) and man pages.

## mclust 5.3 (2017-05)

CRAN release: 2017-05-21

- Added
  [`gmmhd()`](https://mclust-org.github.io/mclust/reference/gmmhd.md)
  function and relative methods.
- Added
  [`MclustDRsubsel()`](https://mclust-org.github.io/mclust/reference/MclustDRsubsel.md)
  function and relative methods.
- Added option to use subset in the hierarchical initialization step
  when a noise component is present.
- [`plot.clustCombi()`](https://mclust-org.github.io/mclust/reference/plot.clustCombi.md)
  presents a menu in interactive sessions, no more need of data for
  classification plots but extract the data from the `clustCombi`
  object.
- Added
  [`combiTree()`](https://mclust-org.github.io/mclust/reference/combiTree.md)
  plot for `clustCombi` objects.
- [`clPairs()`](https://mclust-org.github.io/mclust/reference/clPairs.md)
  now produces a single scatterplot in the bivariate case.
- Fixed a bug in
  [`imputeData()`](https://mclust-org.github.io/mclust/reference/imputeData.md)
  when seed is provided. Now if a seed is provided the data matrix is
  reproducible.
- in
  [`imputeData()`](https://mclust-org.github.io/mclust/reference/imputeData.md)
  and
  [`imputePairs()`](https://mclust-org.github.io/mclust/reference/imputePairs.md)
  some name of arguments have been modified to be coherent with the rest
  of the package.
- Added functions
  [`matchCluster()`](https://mclust-org.github.io/mclust/reference/imputeData.md)
  and
  [`majorityVote()`](https://mclust-org.github.io/mclust/reference/majorityVote.md).
- Rewrite of print and summary methods for `clustCombi` class objects.
- Added
  [`clustCombiOptim()`](https://mclust-org.github.io/mclust/reference/clustCombiOptim.md).
- Fixed a bug in
  [`randomPairs()`](https://mclust-org.github.io/mclust/reference/hcRandomPairs.md)
  when nrow of input data is odd.
- Fixed a bug in
  [`plotDensityMclust2()`](https://mclust-org.github.io/mclust/reference/plot.densityMclust.md),
  [`plotDensityMclustd()`](https://mclust-org.github.io/mclust/reference/plot.densityMclust.md)
  and
  [`surfacePlot()`](https://mclust-org.github.io/mclust/reference/surfacePlot.md)
  when a noise component is present.

## mclust 5.2.3 (2017-03)

CRAN release: 2017-03-13

- Added native routine registration for Fortran code.
- Fixed lowercase argument PACKAGE in
  [`.Fortran()`](https://rdrr.io/r/base/Foreign.html) calls.

## mclust 5.2.2 (2017-01)

CRAN release: 2017-01-22

- Fixed a bug in rare case when performing an extra M step at the end of
  EM algorithm.

## mclust 5.2.1 (2017-01)

CRAN release: 2017-01-03

- Replaced `structure(NULL, *)` with `structure(list(), *)`

## mclust 5.2 (2016-03)

CRAN release: 2016-03-31

- Added argument `x` to
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md)
  to use BIC values from previous computations to avoid recomputing for
  the same models. The same argument and functionality was already
  available in
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).
- Added argument `x` to
  [`mclustICL()`](https://mclust-org.github.io/mclust/reference/mclustICL.md)
  to use ICL values from previous computations to avoid recomputing for
  the same models.
- Fixed a bug on
  [`plot.MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/plot.MclustBoostrap.md)
  for the `"mean"` and `"var"` in the univariate case.
- Fixed uncertainty plots.
- Added functions
  [`as.Mclust()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md)
  and
  [`as.densityMclust()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md)
  to convert object to specific mclust classes.
- Solved a numerical accuracy problem in
  [`qclass()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md)
  when the scale of `x` is (very) large by making the tolerance eps
  scale dependent.
- Use transpose subroutine instead of non-Fortran 77 TRANSPOSE function
  in `mclustaddson.f`.
- Fixed
  [`predict.Mclust()`](https://mclust-org.github.io/mclust/reference/predict.Mclust.md)
  and
  [`predict.MclustDR()`](https://mclust-org.github.io/mclust/reference/predict.MclustDR.md)
  by implementing a more efficient and accurate algorithm for computing
  the densities.

## mclust 5.1 (2015-10)

CRAN release: 2015-10-27

- Fixed slow convergence for VVE and EVE models.
- Fixed a bug in orientation for model VEE.
- Added an extra M-step and parameters update in
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md)
  call via
  [`summaryMclustBIC()`](https://mclust-org.github.io/mclust/reference/summary.mclustBIC.md).

## mclust 5.0.2 (2015-07)

CRAN release: 2015-07-08

- Added option to
  [`MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md)
  for using weighted likelihood bootstrap.
- Added a plot method for `MclustBootstrap` objects.
- Added
  [`errorBars()`](https://mclust-org.github.io/mclust/reference/errorBars.md)
  function.
- Added
  [`clPairsLegend()`](https://mclust-org.github.io/mclust/reference/clPairs.md)
  function.
- Added
  [`covw()`](https://mclust-org.github.io/mclust/reference/covw.md)
  function.
- Fixed rescaling of mixing probabilities in new models.
- Bug fixes.

## mclust 5.0.1 (2015-04)

CRAN release: 2015-04-22

- Fixed bugs.
- Added print method for `hc` objects.

## mclust 5.0.0 (2015-03)

CRAN release: 2015-04-09

- Added the four missing models (EVV, VEE, EVE, VVE) to the mclust
  family. A noise component is allowed, but no prior is available.
- Added
  [`mclustBootstrapLRT()`](https://mclust-org.github.io/mclust/reference/mclustBootstrapLRT.md)
  function (and corresponding print and plot methods) for selecting the
  number of mixture components based on the sequential bootstrap
  likelihood ratio test.
- Added
  [`MclustBootstrap()`](https://mclust-org.github.io/mclust/reference/MclustBootstrap.md)
  function (and corresponding print and summary methods) for performing
  bootstrap inference. This provides standard errors for parameters and
  confidence intervals.
- Added `"A quick tour of mclust"` vignette as html generated using
  rmarkdown and knitr. Older vignettes are included as other
  documentation for the package.
- Modified arguments to
  [`mvn2plot()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md)
  to control colour, lty, lwd, and pch of ellipses and mean point.
- Added functions
  [`emX()`](https://mclust-org.github.io/mclust/reference/emE.md),
  [`emXII()`](https://mclust-org.github.io/mclust/reference/emE.md),
  [`emXXI()`](https://mclust-org.github.io/mclust/reference/emE.md),
  [`emXXX()`](https://mclust-org.github.io/mclust/reference/emE.md),
  [`cdensX()`](https://mclust-org.github.io/mclust/reference/cdensE.md),
  [`cdensXII()`](https://mclust-org.github.io/mclust/reference/cdensE.md),
  [`cdensXXI()`](https://mclust-org.github.io/mclust/reference/cdensE.md),
  and
  [`cdensXXX()`](https://mclust-org.github.io/mclust/reference/cdensE.md),
  to deal with single-component cases, so calling the em function works
  even if `G = 1`.
- Small changes to
  [`icl()`](https://mclust-org.github.io/mclust/reference/icl.md), now
  it is a generic method, with specialized methods for `Mclust` and
  `MclustDA` objects.
- Fixed bug for transformations in the initialization step when some
  variables are constant (i.e. the variance is zero) or a
  one-dimensional data is provided.
- Changed the order of arguments in
  [`hc()`](https://mclust-org.github.io/mclust/reference/hc.md) (and all
  the functions calling it).
- Small modification to `CITATION` file upon request of CRAN
  maintainers.
- Various bug fixes.

## mclust 4.4 (2014-09)

CRAN release: 2014-09-16

- Added option for using transformation of variables in the hierarchical
  initialization step.
- Added
  [`quantileMclust()`](https://mclust-org.github.io/mclust/reference/cdfMclust.md)
  for computing the quantiles from a univariate Gaussian mixture
  distribution.
- Fixed bugs on
  [`summaryMclustBIC()`](https://mclust-org.github.io/mclust/reference/summary.mclustBIC.md),
  [`summaryMclustBICn()`](https://mclust-org.github.io/mclust/reference/summary.mclustBIC.md),
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md)
  to return a matrix of 1s on a single column for `z` even in the case
  of `G = 1`. This is to avoid error on some plots.
- Moved pdf files (previously included as vignettes) to `inst/doc` with
  corresponding `index.html`.

## mclust 4.3 (2014-03)

CRAN release: 2014-03-31

- Fixed bug for
  [`logLik.MclustDA()`](https://mclust-org.github.io/mclust/reference/logLik.MclustDA.md)
  in the univariate case.
- Added argument `"what"` to
  [`predict.densityMclust()`](https://mclust-org.github.io/mclust/reference/predict.densityMclust.md)
  function for choosing what to retrieve, the mixture density or
  component density.
- [`hc()`](https://mclust-org.github.io/mclust/reference/hc.md) function
  has an additional parameter to control if the original variables or a
  transformation of them should be used for hierarchical clustering.
- Added `"hcUse"` argument in
  [`mclust.options()`](https://mclust-org.github.io/mclust/reference/mclust.options.md)
  to be passed as default to
  [`hc()`](https://mclust-org.github.io/mclust/reference/hc.md).
- Added the storing of original data (and class for classification
  models) in the object returned by the main functions.
- Added component `hypvol` to `Mclust` object which provide the
  hypervolume of the noise component when required, otherwise is set to
  `NA`.
- Added a warning when prior is used and BIC returns NAs.
- Fixed bugs in
  [`summary.Mclust()`](https://mclust-org.github.io/mclust/reference/summary.Mclust.md),
  [`print.summary.Mclust()`](https://mclust-org.github.io/mclust/reference/summary.Mclust.md),
  [`plot.Mclust()`](https://mclust-org.github.io/mclust/reference/plot.Mclust.md)
  and [`icl()`](https://mclust-org.github.io/mclust/reference/icl.md) in
  the case of presence of a noise component.
- Fixed bug on some plots in
  [`plot.MclustDR()`](https://mclust-org.github.io/mclust/reference/plot.MclustDR.md)
  which requires [`plot.new()`](https://rdrr.io/r/graphics/frame.html)
  before calling
  [`plot.window()`](https://rdrr.io/r/graphics/plot.window.html).
- Fixed a bug in
  [`MclustDR()`](https://mclust-org.github.io/mclust/reference/MclustDR.md)
  for the one-dimensional case.
- Corrections to `Mclust` man page.
- Various small bug fixes.

## mclust 4.2 (2013-07)

CRAN release: 2013-07-19

- Fixed bug in `sim*()` functions when no obs are assigned to a
  component.
- [`MclustDA()`](https://mclust-org.github.io/mclust/reference/MclustDA.md)
  allows to fit a single class model.
- Fixex bug in
  [`summary.Mclust()`](https://mclust-org.github.io/mclust/reference/summary.Mclust.md)
  when a subset is used for initialization.
- Fixed a bug in the function
  [`qclass()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md)
  when ties are present in quantiles, so it always return the required
  number of classes.
- Various small bug fixes.

## mclust 4.1 (2013-04)

CRAN release: 2013-05-01

- Added [`icl()`](https://mclust-org.github.io/mclust/reference/icl.md)
  function for computing the integrated complete-data likelihood.
- Added
  [`mclustICL()`](https://mclust-org.github.io/mclust/reference/mclustICL.md)
  function with associated print and plot methods.
- [`print.mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)
  shows also the top models based on BIC.
- Modified
  [`summary.Mclust()`](https://mclust-org.github.io/mclust/reference/summary.Mclust.md)
  to return also the icl.
- Rewrite of
  [`adjustedRandIndex()`](https://mclust-org.github.io/mclust/reference/adjustedRandIndex.md)
  function. This version is more efficient for large vectors.
- Updated help for
  [`adjustedRandIndex()`](https://mclust-org.github.io/mclust/reference/adjustedRandIndex.md).
- Modifications to
  [`MclustDR()`](https://mclust-org.github.io/mclust/reference/MclustDR.md)
  and its summary method.
- Changed behavior of `plot.MclustDR(..., what = "contour")`.
- Improved plot of uncertainty for
  `plot.MclustDR(..., what = "boundaries")`.
- Corrected a bug for malformed GvHD data.
- Corrected version of
  [`qclass()`](https://mclust-org.github.io/mclust/reference/mclust-internal.md)
  for selecting initial values in case of 1D data when successive
  quantiles coincide.
- Corrected version of plot BIC values when only a single G-component
  models are fitted.
- Various bug fixes.

## mclust 4.0 (2012-08)

CRAN release: 2012-08-09

- Added new summary and print methods for
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md).
- Added new summary and print methods for
  [`densityMclust()`](https://mclust-org.github.io/mclust/reference/densityMclust.md).
- Included
  [`MclustDA()`](https://mclust-org.github.io/mclust/reference/MclustDA.md)
  function and methods.
- Included
  [`MclustDR()`](https://mclust-org.github.io/mclust/reference/MclustDR.md)
  function and methods.
- Included
  [`me.weighted()`](https://mclust-org.github.io/mclust/reference/me.weighted.md)
  function.
- Restored hierarchical clustering capability for the EEE model (hcEEE).
- Included vignettes for mclust version 4 from Technical Report No.  597
  and for using weights in mclust.
- Adoption of GPL (\>= 2) license.

## mclust 3.5 (2012-07)

CRAN release: 2012-07-22

- Added
  [`summary.Mclust()`](https://mclust-org.github.io/mclust/reference/summary.Mclust.md).
- New functions for plotting and summarizing density estimation.
- Various bug fixes.
- Added
  [`clustCombi()`](https://mclust-org.github.io/mclust/reference/clustCombi.md)
  and related functions (code and doc provided by Jean-Patrick Baudry).
- Bug fix: variable names lost when G = 1.

## mclust 3.4.11 (2012-01)

CRAN release: 2012-01-07

- Added `NAMESPACE`.

## mclust 3.4.10 (2011-05)

CRAN release: 2011-05-30

- Removed intrinsic gamma-

## mclust 3.4.9 (2011-05)

CRAN release: 2011-05-28

- Fixed
  [`hypvol()`](https://mclust-org.github.io/mclust/reference/hypvol.md)
  function to avoid overflow.
- Fixed
  [`hypvol()`](https://mclust-org.github.io/mclust/reference/hypvol.md)
  help file value description.
- Removed unused variables and tabs from source code.
- Switched to intrinsic gamma in source code.
- Fixed default warning in estepVEV and mstepVEV.

## mclust 3.4.8 (2010-12)

CRAN release: 2010-12-12

- Fixed output when G = 1 (it had NA for the missing `z` component).

## mclust 3.4.7 (2010-10)

CRAN release: 2010-10-23

- Removed hierarchical clustering capability for the `EEE` model
  (hcEEE).
- The R 2.12.0 build failed due to a 32-bit Windows compiler error,
  forcing removal of the underlying Fortran code for hcEEE from the
  package, which does not contain errors and compiles on other
  platforms.

## mclust 3.4.6 (2010-08)

CRAN release: 2010-08-10

- Added description of parameters output component to `Mclust` and
  `summary.mclustBIC` help files.

## mclust 3.4.5 (2010-07)

CRAN release: 2010-07-24

- Added
  [`densityMclust()`](https://mclust-org.github.io/mclust/reference/densityMclust.md)
  function.

## mclust 3.4.4 (2010-04)

CRAN release: 2010-04-08

- Fixed bug in covariance matrix output for EEV and VEV models.

## mclust 3.4.3 (2010-02)

CRAN release: 2010-02-20

- Bug fixes.

## mclust 3.4.2 (2010-02)

CRAN release: 2010-02-13

- Moved CITATION to inst and used standard format
- BibTex entries are in inst/cite.
- Fixed bug in handling missing classes in
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md).
- Clarified license wording.

## mclust 3.4.1 (2010-01)

CRAN release: 2010-01-22

- Corrected output description in `mclustModel` help file.
- Updated mclust manual reference to show revision.

## mclust 3.4 (2009-12)

CRAN release: 2009-12-16

- Updated `defaultPrior` help file.
- Added utility functions for imputing missing data with the mix
  package.
- Changed default max to number of mixture components in each class from
  9 to 3.

## mclust 3.3.2 (2009-10)

CRAN release: 2009-10-13

- Fixed problems with in `mclustOptions` help file

## mclust 3.3.1 (2009-06)

CRAN release: 2009-07-01

- Fixed
  [`plot.mclustBIC()`](https://mclust-org.github.io/mclust/reference/plot.mclustBIC.md)
  and
  [`plot.Mclust()`](https://mclust-org.github.io/mclust/reference/plot.Mclust.md)
  to handle `modelNames`.
- Changed “orientation” for VEV, VVV models to be consistent with R
  [`eigen()`](https://rdrr.io/r/base/eigen.html) and the literature
- Fixed some problems including doc for the noise option.
- Updated the
  [`unmap()`](https://mclust-org.github.io/mclust/reference/unmap.md)
  function to optionally include missing groups.

## mclust 3.3 (2009-06)

- Fixed bug in the `"errors"` option for
  [`randProj()`](https://mclust-org.github.io/mclust/reference/randProj.md).
- Fixed boundary cases for the `"noise"` option.

## mclust 3.2 (2009-04)

CRAN release: 2009-04-29

- Added permission for CRAN distribution to LICENSE.
- Fixed problems with help files found by new parser.
- Changed PKG_LIBS order in src/Makevars.
- Fixed
  [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md)
  to handle sampling in data expression in call.

## mclust 3.1.10 (2008-11)

- Added `EXPR = to` all switch functions that didn’t already have it.

## mclust 3.1.9 (2008-10)

- Added `pro` component to parameters in
  [`dens()`](https://mclust-org.github.io/mclust/reference/dens.md) help
  file.
- Fixed some problems with the noise option.

## mclust 3.1.1 (2007-03)

- Default seed changed in `sim*()` functions.
- Added model name check to various functions.
- Otherwise backward compatible with version 3.0

## mclust 3.1 (2007-01)

- Most plotting functions changed to use color.
- [`Mclust()`](https://mclust-org.github.io/mclust/reference/Mclust.md)
  and
  [`mclustBIC()`](https://mclust-org.github.io/mclust/reference/mclustBIC.md)
  fixed to work with G=1
- Otherwise backward compatible with version 3.0.

## mclust 3.0 (2006-10)

- New functionality added, including conjugate priors for Bayesian
  regularization.
- Backward compatibility is not guaranteed since the implementation of
  some functions has changed to make them easier to use or maintain.
