# Template for variance specification for parameterized Gaussian mixture models

Specification of variance parameters for the various types of Gaussian
mixture models.

## Usage

``` r
mclustVariance(modelName, d = NULL, G = 2)
```

## Arguments

- modelName:

  A character string specifying the model.

- d:

  A integer specifying the dimension of the data.

- G:

  An integer specifying the number of components in the mixture model.

## Details

The `variance` component in the `parameters` list from the output to
e.g. `me` or `mstep` or input to e.g. `estep` may contain one or more of
the following arguments, depending on the model:

- `modelName`:

  A character string indicating the model.

- `d`:

  The dimension of the data.

- `G`:

  The number of components in the mixture model.

- `sigmasq`:

  for the one-dimensional models (`"E"`, `"V"`) and spherical models
  (`"EII"`, `"VII"`). This is either a vector whose *k*th component is
  the variance for the *k*th component in the mixture model (`"V"` and
  `"VII"`), or a scalar giving the common variance for all components in
  the mixture model (`"E"` and `"EII"`).

- `Sigma`:

  For the equal variance models `"EII"`, `"EEI"`, and `"EEE"`. A *d* by
  *d* matrix giving the common covariance for all components of the
  mixture model.

- `cholSigma`:

  For the equal variance model `"EEE"`. A *d* by *d* upper triangular
  matrix giving the Cholesky factor of the common covariance for all
  components of the mixture model.

- `sigma`:

  For all multidimensional mixture models. A *d* by *d* by *G* matrix
  array whose `[,,k]`th entry is the covariance matrix for the *k*th
  component of the mixture model.

- `cholsigma`:

  For the unconstrained covariance mixture model `"VVV"`. A *d* by *d*
  by *G* matrix array whose `[,,k]`th entry is the upper triangular
  Cholesky factor of the covariance matrix for the *k*th component of
  the mixture model.

- `scale`:

  For diagonal models `"EEI"`, `"EVI"`, `"VEI"`, `"VVI"` and
  constant-shape models `"EEV"` and `"VEV"`. Either a *G*-vector giving
  the scale of the covariance (the *d*th root of its determinant) for
  each component in the mixture model, or a single numeric value if the
  scale is the same for each component.

- `shape`:

  For diagonal models `"EEI"`, `"EVI"`, `"VEI"`, `"VVI"` and
  constant-shape models `"EEV"` and `"VEV"`. Either a *G* by *d* matrix
  in which the *k*th column is the shape of the covariance matrix
  (normalized to have determinant 1) for the *k*th component, or a
  *d*-vector giving a common shape for all components.

- `orientation`:

  For the constant-shape models `"EEV"` and `"VEV"`. Either a *d* by *d*
  by *G* array whose `[,,k]`th entry is the orthonomal matrix whose
  columns are the eigenvectors of the covariance matrix of the *k*th
  component, or a *d* by *d* orthonormal matrix if the mixture
  components have a common orientation. The `orientation` component is
  not needed in spherical and diagonal models, since the principal
  components are parallel to the coordinate axes so that the orientation
  matrix is the identity.

In all cases, the value `-1` is used as a placeholder for unknown
nonzero entries.
