# Structure search and estimation for Gaussian graphical models

Graph structure search and estimation for Gaussian covariance and
concentration graph models.

## Usage

``` r
searchGGM(data = NULL,
          S = NULL, n = NULL,
          model = c("covariance", "concentration"),
          search = c("step-forw", "step-back", "ga"),
          penalty = c("bic", "ebic", "erdos", "power"),
          beta = NULL,
          start = NULL,
          regularize = FALSE, regHyperPar = NULL,
          ctrlSTEP = controlSTEP(), ctrlGA = controlGA(), ctrlICF = controlICF(),
          parallel = FALSE,
          verbose = FALSE, ...)
```

## Arguments

- data:

  A dataframe or matrix, where rows correspond to observations and
  columns to variables. Categorical variables are not allowed.

- S:

  The sample covariance matrix of the data. If `S = NULL`, the maximum
  likelihood estimate of the covariance matrix is used in the estimation
  of the graphical model.

- n:

  The number of observations. If `data = NULL` and `S` is provided in
  input, `n` must be provided in input as well.

- model:

  The type of Gaussian graphical model. Default is `"covariance"`. See
  "Details".

- search:

  The type of structure search algorithm. If `search = "step-forw"`, a
  greedy forward-stepwise search is used to find the optimal graph
  association structure. If `search = "step-back"`, a greedy
  backward-stepwise search is implemented. If `search = "ga"` a
  stochastic search based on a genetic algorithm is employed. Default is
  `"step-forw"`.

- penalty:

  The penalty function used to define a criterion for scoring the
  candidate graph configurations. Default is `"bic"`. See "Details" and
  [`penalty`](penalty.md).

- beta:

  The hyperparameter of the penalty function. See "Details" and
  [`penalty`](penalty.md).

- start:

  A starting matrix for the estimation algorithm. If `NULL`, the
  starting value is the diagonal sample covariance matrix. Used only
  when `model = "covariance"`.

- regularize:

  A logical argument indicating if Bayesian regularization should be
  performed. Default to `FALSE`. Used only when `model = "covariance"`.

- regHyperPar:

  A list of hyper parameters for Bayesian regularization. Only used when
  `regularization = TRUE`; see also
  [`controlREG`](control-parameters.md).

- ctrlSTEP:

  A list of control parameters used in the stepwise search; see also
  [`controlSTEP`](control-parameters.md).

- ctrlGA:

  A list of control parameters for the genetic algorithm; see also
  [`controlGA`](control-parameters.md).

- ctrlICF:

  A list of control parameters employed in the algorithm for estimation
  of graphical model parameters; see also
  [`controlICF`](control-parameters.md).

- parallel:

  A logical argument indicating if parallel computation should be used
  for structure search. If TRUE, all the available cores are used. The
  argument could also be set to a numeric integer value specifying the
  number of cores to be employed.

- verbose:

  A logical argument controlling whether iterations of the structure
  searching and estimation procedure need to be shown or not.

- ...:

  Additional internal arguments not to be provided by the user.

## Value

An object of class `'fitGGM'` containing the optimal estimated marginal
or conditional independence Gaussian graphical model.

The output is a list containing:

- sigma:

  The estimated covariance matrix.

- omega:

  The estimated concentration (inverse covariance) matrix.

- graph:

  The adjacency matrix corresponding to the optimal marginal or
  conditional independence graph.

- model:

  Estimated model type, whether `"covariance"` or `"concentration"`.

- loglikPen:

  Value of the maximized penalized log-likelihood.

- loglik:

  Value of the maximized log-likelihood.

- nPar:

  Number of estimated parameters.

- n:

  Number of observations.

- V:

  Number of variables, corresponding to the number of nodes in the
  graph.

- penalty:

  The type of penalty on the graph structure.

- search:

  The search method used for graph structure search.

- GA:

  An object of class `'ga-class'` with information about the genetic
  algorithm. Only present when `search = "ga"`. See
  [`ga`](https://github.com/luca-scr/GA/reference/ga.html).

## Details

The function performs graph association structure search and maximum
penalized likelihood estimation of the optimal Gaussian graphical model
given the data provided in input.

A Gaussian covariance graph model is estimated if
`model = "covariance"`, while estimation of a Gaussian covariance graph
model is performed if `model = "concentration"`. A Gaussian covariance
graph model postulates that some variables are marginally independent
according to the inferred graph structure. On the other hand, in a
Gaussian concentration graph model, variables are conditionally
independent given their neighbors in the inferred graph. See also
[`fitGGM`](fitGGM.md).

Search for the optimal graph structure and parameter estimation is
carried out by maximization of a Gaussian penalized likelihood, given as
follows: \$\$ \mbox{Covariance:}\~\~~ \mbox{argmax}\_{\Sigma,A} \~~
\ell(X \| \Sigma, A) - P(A, \beta) \quad \Sigma \in C^{+}\_G(A)\$\$ \$\$
\mbox{Concentration:}\~\~~ \mbox{argmax}\_{\Omega,A} \~~ \ell(X \|
\Omega, A) - P(A, \beta) \quad \Omega \in C^{+}\_G(A)\$\$ where
\\(C^{+}\_G(A)\\) is the collection of sparse positive definite matrices
whose zero patterns are given by graph \\(G\\) represented by the
adjacency matrix \\(A\\).

The penalty function \\(P(A, \beta)\\) depends on the structure of graph
\\(G\\) through the adjacency matrix \\(A\\) and a parameter
\\(\beta\\); see [`penalty`](penalty.md) on how to specify the
penalization term and for further information.

For this type of penalized log-likelihood, graph structure search and
parameter estimation is a maximization combinatorial problem. For a
given candidate structure (i.e. adjacency matrix), association
parameters in the covariance or concentration matrix are estimated using
the estimation algorithms implemented in [`fitGGM`](fitGGM.md).
Regarding structure search, this can be carried out either using a
greedy forward-stepwise or a greedy backward-stepwise algorithm, by
setting `search = "step-forw"` or `search = "step-back"` respectively.
Alternatively, a stochastic search via genetic algorithm can be used by
setting `search = "ga"`. The procedure for the forward stepwise search
is described in Fop et al. (2018), and the backward is implemented in a
similar way; the genetic algorithm procedure relies on the
[`GA`](https://github.com/luca-scr/GA/reference/GA-package.html)
package. All the structure searching methods can be run in parallel on a
multi-core machine by setting the argument `parallel = TRUE`.

## References

Fop, M., Murphy, T.B., and Scrucca, L. (2019) Model-based clustering
with sparse covariance matrices. \*\*Statistics and Computing\*\*, 29:4,
791–819.

Scrucca, L. (2017). On some extensions to GA package: Hybrid
optimisation, parallelisation and islands evolution. *The R Journal*,
9(1), 187-206.

Scrucca, L. (2013). GA: A package for genetic algorithms in R. *Journal
of Statistical Software*, 53(4), 1-3.

## Examples

``` r

# fit covariance graph model with default forward-stepwise search
data(mtcars)
x <- mtcars[,c(1,3:7)]
mod1 <- searchGGM(x, model = "covariance")
mod1
#> 
#> 'fitGGM' object containing: 
#>  [1] "omega"     "graph"     "model"     "loglikPen" "loglik"    "nPar"     
#>  [7] "n"         "V"         "penalty"   "search"   
plot(mod1)

#
# prefer a sparser model
mod2 <- searchGGM(x, model = "covariance", penalty = "ebic")
mod2
#> 
#> 'fitGGM' object containing: 
#>  [1] "omega"     "graph"     "model"     "loglikPen" "loglik"    "nPar"     
#>  [7] "n"         "V"         "penalty"   "search"   
plot(mod2)



# fit concentration graph model with backward-stepwise structure search
# with a covariance matrix in input
data(ability.cov)
mod3 <- searchGGM(S = ability.cov$cov, n = ability.cov$n.obs,
                  model = "concentration", search = "step-back")
mod3
#> 
#> 'fitGGM' object containing: 
#>  [1] "omega"     "graph"     "model"     "loglikPen" "loglik"    "nPar"     
#>  [7] "n"         "V"         "penalty"   "search"   
mod3$graph
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    0    1    1    0    1    0
#> [2,]    1    0    1    0    0    0
#> [3,]    1    1    0    1    0    0
#> [4,]    0    0    1    0    0    0
#> [5,]    1    0    0    0    0    1
#> [6,]    0    0    0    0    1    0
mod3$omega
#>             general     picture      blocks        maze     reading       vocab
#> general  0.08115994 -0.02581133 -0.01049276  0.00000000 -0.02377784  0.00000000
#> picture -0.02581133  0.23330570 -0.02235247  0.00000000  0.00000000  0.00000000
#> blocks  -0.01049276 -0.02235247  0.01338466 -0.01271862  0.00000000  0.00000000
#> maze     0.00000000  0.00000000 -0.01271862  0.09810766  0.00000000  0.00000000
#> reading -0.02377784  0.00000000  0.00000000  0.00000000  0.06033143 -0.02510099
#> vocab    0.00000000  0.00000000  0.00000000  0.00000000 -0.02510099  0.01977790
plot(mod3)


if (FALSE) { # \dontrun{

# generate data from a Markov model
n <- 1000
V <- 20
dat <- matrix(NA, n, V)
dat[,1] <- rnorm(n)
for ( j in 2:V ) dat[,j] <- dat[,j-1] + rnorm(n, sd = 0.5)
mod4 <- searchGGM(data = dat, model = "concentration")        # recover the model
plot(mod4, what = "adjacency")

} # }
```
