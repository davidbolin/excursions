# Sequential estimation of Gaussian integrals

`gaussint` is used for calculating \\n\\-dimensional Gaussian integrals
\$\$\int_a^b \frac{\|Q\|^{1/2}}{(2\pi)^{n/2}}
\exp(-\frac1{2}(x-\mu)^{T}Q(x-\mu)) dx\$\$ A limit value \\lim\\ can be
used to stop the integration if the sequential estimate goes below the
limit, which can result in substantial computational savings in cases
when one only is interested in testing if the integral is above the
limit value. The integral is calculated sequentially, and estimates for
all subintegrals are also returned.

## Usage

``` r
gaussint(
  mu,
  Q.chol,
  Q,
  a,
  b,
  lim = 0,
  n.iter = 10000,
  ind,
  use.reordering = c("natural", "sparsity", "limits"),
  max.size,
  max.threads = 0,
  seed
)
```

## Arguments

- mu:

  Expectation vector for the Gaussian distribution.

- Q.chol:

  The Cholesky factor of the precision matrix (optional).

- Q:

  Precision matrix for the Gaussian distribution. If Q is supplied but
  not Q.chol, the cholesky factor is computed before integrating.

- a:

  Lower limit in integral.

- b:

  Upper limit in integral.

- lim:

  If this argument is used, the integration is stopped and 0 is returned
  if the estimated value goes below \\lim\\.

- n.iter:

  Number or iterations in the MC sampler that is used for approximating
  probabilities. The default value is 10000.

- ind:

  Indices of the nodes that should be analyzed (optional).

- use.reordering:

  Determines what reordering to use:

  "natural"

  :   No reordering is performed.

  "sparsity"

  :   Reorder for sparsity in the cholesky factor (MMD reordering is
      used).

  "limits"

  :   Reorder by moving all nodes with a=-Inf and b=Inf first and then
      reordering for sparsity (CAMD reordering is used).

- max.size:

  The largest number of sub-integrals to compute. Default is the total
  dimension of the distribution.

- max.threads:

  Decides the number of threads the program can use. Set to 0 for using
  the maximum number of threads allowed by the system (default).

- seed:

  The random seed to use (optional).

## Value

A list with elements

- P :

  Value of the integral.

- E :

  Estimated error of the P estimate.

- Pv :

  A vector with the estimates of all sub-integrals.

- Ev :

  A vector with the estimated errors of the Pv estimates.

## Details

The function uses sequential importance sampling to estimate the
Gaussian integral, and returns all computed sub-integrals. This means
that if, for example, the function is used to compute \\P(x\>0)\\ for an
n-dimensional Gaussian variable \\x\\, then all integrals
\\P(x_1\>0,\ldots,x_i\>0)\\ for \\i=1,\ldots,n\\ are computed.

If one is only interested in whether \\P(x\>0)\>\alpha\\ or not, then
one can stop the integration as soon as
\\P(x_1\>0,\ldots,x_i\>0)\<\alpha\\. This can save a lot of computation
time if \\P(x_1\>0,\ldots,x_i\>0)\< \alpha\\ for \\i\\ much smaller than
\\n\\. This limit value is specified by the `lim` argument.

Which reordering to use depends on what the purpose of the calculation
is and what the integration limits are. However, in general the `limits`
reordering is typically most appropriate since this combines sparisty
(which improves accuracy and reduces computational cost) with automatic
handling of dimensions with limits `a=-Inf` and `b=Inf`, which do not
affect the probability but affect the computation time if they are not
handled separately.

## References

Bolin, D. and Lindgren, F. (2015) *Excursion and contour uncertainty
regions for latent Gaussian models*, JRSS-series B, vol 77, no 1, pp
85-106.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, vol 86, no 1, pp 1-20.

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
## Create mean and a tridiagonal precision matrix
n <- 11
mu.x <- seq(-5, 5, length = n)
Q.x <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
## Calculate the probability that the variable is between mu-3 and mu+3
prob <- gaussint(mu = mu.x, Q = Q.x, a = mu.x - 3, b = mu.x + 3, max.threads = 2)
prob$P
#> [1] 0.9680017
```
