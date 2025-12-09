# Excursion sets and contour credible regions for latent Gaussian models

Excursion sets and contour credible regions for latent Gaussian models
calculated using the INLA method.

## Usage

``` r
excursions.inla(
  result.inla,
  stack,
  name = NULL,
  tag = NULL,
  ind = NULL,
  method,
  alpha = 1,
  F.limit,
  u,
  u.link = FALSE,
  type,
  n.iter = 10000,
  verbose = 0,
  max.threads = 0,
  compressed = TRUE,
  seed = NULL,
  prune.ind = FALSE
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

- method:

  Method for handeling the latent Gaussian structure:

  'EB'

  :   Empirical Bayes

  'QC'

  :   Quantile correction

  'NI'

  :   Numerical integration

  'NIQC'

  :   Numerical integration with quantile correction

  'iNIQC'

  :   Improved integration with quantile correction

- alpha:

  Error probability for the excursion set of interest. The default value
  is 1.

- F.limit:

  Error probability for when to stop the calculation of the excursion
  function. The default value is `alpha`, and the value cannot be
  smaller than `alpha`. A smaller value of `F.limit` results in a
  smaller computation time.

- u:

  Excursion or contour level.

- u.link:

  If u.link is TRUE, `u` is assumed to be in the scale of the data and
  is then transformed to the scale of the linear predictor (default
  FALSE).

- type:

  Type of region:

  '\>'

  :   positive excursions

  '\<'

  :   negative excursions

  '!='

  :   contour avoiding function

  '='

  :   contour credibility function

- n.iter:

  Number or iterations in the MC sampler that is used for approximating
  probabilities. The default value is 10000.

- verbose:

  Set to TRUE for verbose mode (optional).

- max.threads:

  Decides the number of threads the program can use. Set to 0 for using
  the maximum number of threads allowed by the system (default).

- compressed:

  If INLA is run in compressed mode and a part of the linear predictor
  is to be used, then only add the relevant part. Otherwise the entire
  linear predictor is added internally (default TRUE).

- seed:

  Random seed (optional).

- prune.ind:

  If `TRUE` and `ind` is supplied, then the result object is pruned to
  contain only the active nodes specified by `ind`.

## Value

`excursions.inla` returns an object of class "excurobj" with the
following elements

- E :

  Excursion set, contour credible region, or contour avoiding set

- F :

  The excursion function corresponding to the set `E` calculated for
  values up to `F.limit`

- G :

  Contour map set. \\G=1\\ for all nodes where the \\mu \> u\\.

- M :

  Contour avoiding set. \\M=-1\\ for all non-significant nodes. \\M=0\\
  for nodes where the process is significantly below `u` and \\M=1\\ for
  all nodes where the field is significantly above `u`. Which values
  that should be present depends on what type of set that is calculated.

- rho :

  Marginal excursion probabilities

- mean :

  Posterior mean

- vars :

  Marginal variances

- meta :

  A list containing various information about the calculation.

## Details

The different methods for handling the latent Gaussian structure are
listed in order of accuracy and computational cost. The `EB` method is
the simplest and is based on a Gaussian approximation of the posterior
of the quantity of interest. The `QC` method uses the same Gaussian
approximation but improves the accuracy by modifying the limits in the
integrals that are computed in order to find the region. The other three
methods are intended for Bayesian models where the posterior
distribution for the quantity of interest is obtained by integrating
over the parameters in the model. The `NI` method approximates this
integration in the same way as is done in INLA, and the `NIQC` and
`iNIQC` methods combine this apprximation with the QC method for
improved accuracy.

If the main purpose of the analysis is to construct excursion or contour
sets for low values of `alpha`, we recommend using `QC` for problems
with Gaussian likelihoods and `NIQC` for problems with non-Gaussian
likelihoods. The reason for this is that the more accurate methods also
have higher computational costs.

## Note

This function requires the `INLA` package, which is not a CRAN package.
See <https://www.r-inla.org/download-install> for easy installation
instructions.

## References

Bolin, D. and Lindgren, F. (2015) *Excursion and contour uncertainty
regions for latent Gaussian models*, JRSS-series B, vol 77, no 1, pp
85-106.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, vol 86, no 1, pp 1-20.

## See also

[`excursions()`](https://davidbolin.github.io/excursions/reference/excursions.md),
[`excursions.mc()`](https://davidbolin.github.io/excursions/reference/excursions.mc.md)

## Author

David Bolin <davidbolin@gmail.com> and Finn Lindgren
<finn.lindgren@gmail.com>

## Examples

``` r
## In this example, we calculate the excursion function
## for a partially observed AR process.
if (FALSE) { # \dontrun{
if (require.nowarnings("INLA")) {
  ## Sample the process:
  rho <- 0.9
  tau <- 15
  tau.e <- 1
  n <- 100
  x <- 1:n
  mu <- 10 * ((x < n / 2) * (x - n / 2) + (x >= n / 2) * (n / 2 - x) + n / 4) / n
  Q <- tau * sparseMatrix(
    i = c(1:n, 2:n), j = c(1:n, 1:(n - 1)),
    x = c(1, rep(1 + rho^2, n - 2), 1, rep(-rho, n - 1)),
    dims = c(n, n), symmetric = TRUE
  )
  X <- mu + solve(chol(Q), rnorm(n))

  ## measure the sampled process at n.obs random locations
  ## under Gaussian measurement noise.
  n.obs <- 50
  obs.loc <- sample(1:n, n.obs)
  A <- sparseMatrix(
    i = 1:n.obs, j = obs.loc, x = rep(1, n.obs),
    dims = c(n.obs, n)
  )
  Y <- as.vector(A %*% X + rnorm(n.obs) / sqrt(tau.e))

  ## Estimate the parameters using INLA
  ef <- list(c(list(ar = x), list(cov = mu)))
  s.obs <- inla.stack(data = list(y = Y), A = list(A), effects = ef, tag = "obs")
  s.pre <- inla.stack(data = list(y = NA), A = list(1), effects = ef, tag = "pred")
  stack <- inla.stack(s.obs, s.pre)
  formula <- y ~ -1 + cov + f(ar, model = "ar1")
  result <- inla(
    formula = formula, family = "normal", data = inla.stack.data(stack),
    control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
    control.compute = list(
      config = TRUE,
      return.marginals.predictor = TRUE
    )
  )

  ## calculate the level 0 positive excursion function
  res.qc <- excursions.inla(result,
    stack = stack, tag = "pred", alpha = 0.99, u = 0,
    method = "QC", type = ">", max.threads = 2
  )
  ## plot the excursion function and marginal probabilities
  plot(res.qc$rho,
    type = "l",
    main = "marginal probabilities (black) and excursion function (red)"
  )
  lines(res.qc$F, col = 2)
}
} # }
```
