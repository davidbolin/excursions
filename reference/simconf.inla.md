# Simultaneous confidence regions for latent Gaussian models

`simconf.inla` is used for calculating simultaneous confidence regions
for latent Gaussian models estimated using `INLA`.

## Usage

``` r
simconf.inla(
  result.inla,
  stack,
  name = NULL,
  tag = NULL,
  ind = NULL,
  alpha,
  method = "NI",
  n.iter = 10000,
  verbose = FALSE,
  link = FALSE,
  max.threads = 0,
  compressed = TRUE,
  seed = NULL,
  inla.sample = TRUE
)
```

## Arguments

- result.inla:

  Result object from INLA call.

- stack:

  The stack object used in the INLA call.

- name:

  The name of the component for which to do the calculation. This
  argument should only be used if a stack object is not provided, use
  the tag argument otherwise.

- tag:

  The tag of the component in the stack for which to do the calculation.
  This argument should only be used if a stack object is provided, use
  the name argument otherwise.

- ind:

  If only a part of a component should be used in the calculations, this
  argument specifies the indices for that part.

- alpha:

  Error probability for the region.

- method:

  Method for handling the latent Gaussian structure:

  'EB'

  :   Empirical Bayes (Gaussian approximation of posterior).

  'NI'

  :   Numerical integration (Calculation based on the Gaussian mixture
      approximation of the posterior, as calculated by INLA).

- n.iter:

  Number or iterations in the MC sampler that is used for approximating
  probabilities. The default value is 10000.

- verbose:

  Set to TRUE for verbose mode (optional).

- link:

  Transform output to the scale of the data using the link function as
  defined in the model estimated with INLA (default FALSE).

- max.threads:

  Decides the number of threads the program can use. Set to 0 for using
  the maximum number of threads allowed by the system (default).

- compressed:

  If INLA is run in compressed mode and a part of the linear predictor
  is to be used, then only add the relevant part. Otherwise the entire
  linear predictor is added internally (default TRUE).

- seed:

  Random seed (optional).

- inla.sample:

  Set to TRUE if inla.posterior.sample should be used for the MC
  integration.

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

See
[`simconf()`](https://davidbolin.github.io/excursions/reference/simconf.md)
for details.

## Note

This function requires the `INLA` package, which is not a CRAN package.
See <https://www.r-inla.org/download-install> for easy installation
instructions.

## References

Bolin et al. (2015) *Statistical prediction of global sea level from
global temperature*, Statistica Sinica, vol 25, pp 351-367.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, vol 86, no 1, pp 1-20.

## See also

[`simconf()`](https://davidbolin.github.io/excursions/reference/simconf.md),
[`simconf.mc()`](https://davidbolin.github.io/excursions/reference/simconf.mc.md),
[`simconf.mixture()`](https://davidbolin.github.io/excursions/reference/simconf.mixture.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
if (require.nowarnings("INLA")) {
  n <- 10
  x <- seq(0, 6, length.out = n)
  y <- sin(x) + rnorm(n)
  mu <- 1:n
  result <- inla(y ~ 1 + f(mu, model = "rw2"),
    data = list(y = y, mu = mu), verbose = FALSE,
    control.compute = list(
      config = TRUE,
      return.marginals.predictor = TRUE
    ),
    num.threads = "1:1"
  )

  res <- simconf.inla(
    result,
    name = "mu", alpha = 0.05,
    max.threads = 1, num.threads = "1:1"
  )

  plot(result$summary.random$mu$mean, ylim = c(-2, 2))
  lines(res$a)
  lines(res$b)
  lines(res$a.marginal, col = "2")
  lines(res$b.marginal, col = "2")
}
} # }
```
