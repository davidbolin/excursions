# Simultaneous confidence regions for Gaussian models

`simconf` is used for calculating simultaneous confidence regions for
Gaussian models \\x\\. The function returns upper and lower bounds \\a\\
and \\b\\ such that \\P(a\<x\<b) = 1-\alpha\\.

## Usage

``` r
simconf(
  alpha,
  mu,
  Q,
  n.iter = 10000,
  Q.chol,
  vars,
  ind = NULL,
  verbose = 0,
  max.threads = 0,
  seed = NULL,
  tol = NULL
)
```

## Arguments

- alpha:

  Error probability for the region.

- mu:

  Expectation vector for the Gaussian distribution.

- Q:

  Precision matrix for the Gaussian distribution.

- n.iter:

  Number or iterations in the MC sampler that is used for approximating
  probabilities. The default value is 10000. If `tol` is given, this is
  the maximal number of iterations.

- Q.chol:

  The Cholesky factor of the precision matrix (optional).

- vars:

  Precomputed marginal variances (optional).

- ind:

  Indices of the nodes that should be analyzed (optional).

- verbose:

  Set to TRUE for verbose mode (optional).

- max.threads:

  The number of threads that the program can use. The default, 0, uses
  the default number of threads of OpenMP, which can be set with the
  environment variable `OMP_NUM_THREADS`. The number of threads is at
  most `OMP_THREAD_LIMIT`, and is one if the package was built without
  OpenMP.

- seed:

  Random seed (optional).

- tol:

  Target for the estimated errors of the probabilities that are computed
  to find the band (optional). If `tol` is given, the number of
  iterations is chosen adaptively, using at most `n.iter` iterations,
  see
  [`gaussint()`](https://davidbolin.github.io/excursions/reference/gaussint.md).
  By default, `n.iter` iterations are always used.

## Value

An object of class "excurobj" with elements

- a :

  The lower bound.

- b :

  The upper bound.

- a.marginal :

  The lower bound for pointwise confidence bands.

- b.marginal :

  The upper bound for pointwise confidence bands.

## Details

The pointwise confidence bands are based on the marginal quantiles,
meaning that `a.marignal` is a vector where the ith element equals
\\\mu_i + q\_{\alpha,i}\\ and `b.marginal` is a vector where the ith
element equals \\\mu_i + q\_{1-\alpha,i}\\, where \\\mu_i\\ is the
expected value of the \\x_i\\ and \\q\_{\alpha,i}\\ is the
\\\alpha\\-quantile of \\x_i-\mu_i\\.

The simultaneous confidence band is defined by the lower limit vector
`a` and the upper limit vector `b`, where \\a_i = \mu_i +c q\_{\alpha}\\
and \\b_i = \mu_i + c q\_{1-\alpha}\\, where \\c\\ is a constant
computed such that \\P(a \< x \< b) = 1-\alpha\\.

## References

Bolin et al. (2015) *Statistical prediction of global sea level from
global temperature*, Statistica Sinica, vol 25, pp 351-367.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, vol 86, no 1, pp 1-20.

## See also

[`simconf.inla()`](https://davidbolin.github.io/excursions/reference/simconf.inla.md),
[`simconf.mc()`](https://davidbolin.github.io/excursions/reference/simconf.mc.md),
[`simconf.mixture()`](https://davidbolin.github.io/excursions/reference/simconf.mixture.md)

## Author

David Bolin <davidbolin@gmail.com> and Finn Lindgren
<finn.lindgren@gmail.com>

## Examples

``` r
## Create mean and a tridiagonal precision matrix
n <- 11
mu.x <- seq(-5, 5, length = n)
Q.x <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
## calculate the confidence region
conf <- simconf(0.05, mu.x, Q.x, max.threads = 2)
## Plot the region
plot(mu.x,
  type = "l", ylim = c(-10, 10),
  main = "Mean (black) and confidence region (red)"
)
lines(conf$a, col = 2)
lines(conf$b, col = 2)
```
