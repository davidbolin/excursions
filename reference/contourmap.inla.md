# Contour maps and contour map quality measures for latent Gaussian models

An interface to the `contourmap` function for latent Gaussian models
calculated using the INLA method.

## Usage

``` r
contourmap.inla(
  result.inla,
  stack,
  name = NULL,
  tag = NULL,
  method = "QC",
  n.levels,
  type = c("standard", "pretty", "equalarea"),
  compute = list(F = TRUE, measures = NULL),
  alpha,
  F.limit,
  n.iter = 10000,
  verbose = FALSE,
  max.threads = 0,
  compressed = TRUE,
  seed = NULL,
  ind,
  ...
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

- method:

  Method for handeling the latent Gaussian structure. Currently only
  Empirical Bayes (EB) and Quantile corrections (QC) are supported.

- n.levels:

  Number of levels in contour map.

- type:

  Type of contour map. One of:

  'standard'

  :   Equidistant levels between smallest and largest value of the
      posterior mean (default).

  'pretty'

  :   Equally spaced 'round' values which cover the range of the values
      in the posterior mean.

  'equalarea'

  :   Levels such that different spatial regions are approximately equal
      in size.

- compute:

  A list with quality indices to compute

  'F':

  :   TRUE/FALSE indicating whether the contour map function should be
      computed (default TRUE)

  'measures':

  :   A list with the quality measures to compute ("P0", "P1", "P2") or
      corresponding bounds based only on the marginal probabilities
      ("P0-bound", "P1-bound", "P2-bound")

- alpha:

  Maximal error probability in contour map function (default=1)

- F.limit:

  The limit value for the computation of the F function. F is set to NA
  for all nodes where F\<1-F.limit. Default is F.limit = `alpha`.

- n.iter:

  Number or iterations in the MC sampler that is used for calculating
  the quantities in `compute`. The default value is 10000.

- verbose:

  Set to TRUE for verbose mode (optional)

- max.threads:

  Decides the number of threads the program can use. Set to 0 for using
  the maximum number of threads allowed by the system (default).

- compressed:

  If INLA is run in compressed mode and a part of the linear predictor
  is to be used, then only add the relevant part. Otherwise the entire
  linear predictor is added internally (default TRUE).

- seed:

  Random seed (optional).

- ind:

  If only a part of a component should be used in the calculations, this
  argument specifies the indices for that part (optional).

- ...:

  Additional arguments to the contour map function. See the
  documentation for `contourmap` for details.

## Value

`contourmap.inla` returns an object of class "excurobj" with the same
elements as returned by `contourmap`.

## Details

The INLA approximation of the quantity of interest is in general a
weighted sum of Gaussian distributions with different parameters. If
`method = 'EB'` is used, then the contour map is computed for the mean
of the component in the weighted sum that has parameters with the
highest likelihood. If on the other hand `method='QC'`, then the contour
map is computed for the posterior mean reported by INLA. If the EB
method also is used in INLA, then this reported posterior mean is equal
to the mean of the component with the highest likelihood. Therefore,
`method='EB'` is appropriate if the EB method also is used in INLA, but
`method='QC'` should be used in general.

The `n.levels` contours in the contour map are are placed according to
the argument `type`. A number of quality measures can be computed based
based on the specified contour map and the distribution of the component
of interest. What should be computed is specified using the `compute`
argument. For details on these quanties, see the reference below.

## Note

This function requires the `INLA` package, which is not a CRAN package.
See <https://www.r-inla.org/download-install> for easy installation
instructions.

## References

Bolin, D. and Lindgren, F. (2017) *Quantifying the uncertainty of
contour maps*, Journal of Computational and Graphical Statistics, 26:3,
513-524.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, vol 86, no 1, pp 1-20.

## See also

[`contourmap()`](https://davidbolin.github.io/excursions/reference/contourmap.md),
[`contourmap.mc()`](https://davidbolin.github.io/excursions/reference/contourmap.mc.md),
[`contourmap.colors()`](https://davidbolin.github.io/excursions/reference/contourmap.colors.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
if (require.nowarnings("INLA")) {
  # Generate mesh and SPDE model
  n.lattice <- 10 # increase for more interesting, but slower, examples
  x <- seq(from = 0, to = 10, length.out = n.lattice)
  lattice <- fmesher::fm_lattice_2d(x = x, y = x)
  mesh <- fmesher::fm_rcdt_2d_inla(
    lattice = lattice,
    extend = FALSE,
    refine = FALSE
  )
  spde <- inla.spde2.matern(mesh, alpha = 2)

  # Generate an artificial sample
  sigma2.e <- 0.01
  n.obs <- 100
  obs.loc <- cbind(
    runif(n.obs) * diff(range(x)) + min(x),
    runif(n.obs) * diff(range(x)) + min(x)
  )
  Q <- inla.spde2.precision(spde, theta = c(log(sqrt(0.5)), log(sqrt(1))))
  x <- inla.qsample(Q = Q, num.threads = 1)
  A <- fmesher::fm_basis(mesh = mesh, loc = obs.loc)
  Y <- as.vector(A %*% x + rnorm(n.obs) * sqrt(sigma2.e))

  ## Estimate the parameters using INLA
  mesh.index <- inla.spde.make.index(name = "field", n.spde = spde$n.spde)
  ef <- list(c(mesh.index, list(Intercept = 1)))

  s.obs <- inla.stack(
    data = list(y = Y),
    A = list(A),
    effects = ef,
    tag = "obs"
  )
  s.pre <- inla.stack(
    data = list(y = NA),
    A = list(1),
    effects = ef,
    tag = "pred"
  )
  stack <- inla.stack(s.obs, s.pre)
  formula <- y ~ -1 + Intercept + f(field, model = spde)
  result <- inla(
    formula = formula, family = "normal", data = inla.stack.data(stack),
    control.predictor = list(
      A = inla.stack.A(stack),
      compute = TRUE
    ),
    control.compute = list(
      config = TRUE,
      return.marginals.predictor = TRUE
    ),
    num.threads = 1
  )

  ## Calculate contour map with two levels
  map <- contourmap.inla(result,
    stack = stack, tag = "pred",
    n.levels = 2, alpha = 0.1, F.limit = 0.1,
    max.threads = 1
  )

  ## Plot the results
  cols <- contourmap.colors(map,
    col = heat.colors(100, 1),
    credible.col = grey(0.5, 1)
  )
  image(matrix(map$M[mesh$idx$lattice], n.lattice, n.lattice), col = cols)
}
} # }
```
