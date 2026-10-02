# Connected excursion regions for latent Gaussian models

Connected excursion regions for latent Gaussian models fitted with
`INLA` or `inlabru`. See
[`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md)
for details on the regions.

## Usage

``` r
excursions.regions.inla(
  result.inla,
  stack,
  name = NULL,
  tag = NULL,
  ind = NULL,
  method,
  alpha,
  u,
  u.link = FALSE,
  type,
  graph,
  n.iter = 20000,
  max.regions = Inf,
  min.size = 1,
  growth = c("bound", "rho"),
  n.starts = 10,
  min.prominence = 0,
  verbose = 0,
  max.threads = 0,
  compressed = TRUE,
  seed = NULL,
  prune.ind = FALSE,
  tol = NULL,
  size.tol = 0.001
)
```

## Arguments

- result.inla:

  Result object from an `INLA` or `inlabru` call.

- stack:

  The stack object used in the INLA call.

- name:

  The name of the component for which to do the calculation. This
  argument should only be used if a stack object is not provided, use
  the tag argument otherwise. For `inlabru` results, use
  `name = "APredictor"` together with `ind = bru_index(result, tag)`.

- tag:

  The tag of the component in the stack for which to do the calculation.
  This argument should only be used if a stack object is provided, use
  the name argument otherwise.

- ind:

  If only a part of a component should be used in the calculations, this
  argument specifies the indices for that part.

- method:

  Method for handling the latent Gaussian structure:

  'EB'

  :   Empirical Bayes

  'QC'

  :   Quantile correction

  'NI'

  :   Numerical integration

  'NIQC'

  :   Numerical integration with quantile correction

- alpha:

  Error probability for each region.

- u:

  Excursion level.

- u.link:

  If u.link is TRUE, `u` is assumed to be in the scale of the data and
  is then transformed to the scale of the linear predictor (default
  FALSE).

- type:

  Type of region, `'>'` for positive excursion regions and `'<'` for
  negative excursion regions.

- graph:

  The neighbourhood graph of the nodes, given as a symmetric sparse
  matrix where the non-zero off-diagonal elements are the edges, or as
  an `fm_mesh_2d` object, in which case the vertex graph of the mesh is
  used. The graph can either have one node for each of the nodes
  selected by `ind`, in the same order, or one node for each node of the
  component. The default is the graph of the non-zero elements of the
  precision matrix, which for the linear predictor usually has no edges,
  so the graph should then be provided.

- n.iter:

  Number or iterations in the MC sampler that is used for approximating
  probabilities. The default value is 20000. If `size.tol` or `tol` is
  given, this is the maximal number of iterations.

- max.regions:

  The maximum number of regions to compute.

- min.size:

  The minimum number of nodes of a region.

- growth:

  How the regions are grown, see
  [`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md).

- n.starts:

  The maximum number of start points for growing a region in each
  connected component, see
  [`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md).

- min.prominence:

  The minimum prominence of a local maximum for it to be used as a start
  point, see
  [`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md).

- verbose:

  Set to TRUE for verbose mode (optional).

- max.threads:

  The number of threads that the program can use. The default, 0, uses
  the default number of threads of OpenMP.

- compressed:

  If INLA is run in compressed mode and a part of the linear predictor
  is to be used, then only add the relevant part. Otherwise the entire
  linear predictor is added internally (default TRUE).

- seed:

  Random seed (optional).

- prune.ind:

  If `TRUE` and `ind` is supplied, then the result object is pruned to
  contain only the active nodes specified by `ind`, and the regions are
  given as indices within `ind`.

- tol:

  Target for the estimated errors of the joint probabilities of the
  regions (optional), see
  [`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md).

- size.tol:

  Target for the estimated Monte Carlo error of the size of each region,
  relative to the size, see
  [`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md).

## Value

`excursions.regions.inla` returns a list with the elements

- regions:

  A list with the indices of the regions, largest first. The indices
  refer to the nodes of the component, or to the positions in `ind` if
  `prune.ind = TRUE`.

- P:

  The estimated joint excursion probability of each region.

- P.err:

  The Monte Carlo standard errors of `P`.

- labels:

  A vector with the region of each node, where 0 means that the node is
  not in a region, and `NA` that it is not in `ind`.

- E:

  The largest region, as an indicator vector.

- F:

  The excursion functions of the regions, as a sparse matrix with one
  column for each region, see
  [`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md).

- rho:

  Marginal excursion probabilities.

- mean:

  Posterior mean.

- vars:

  Marginal variances.

- meta:

  A list containing various information about the calculation.

## Details

The methods for handling the latent Gaussian structure are the same as
in
[`excursions.inla()`](https://davidbolin.github.io/excursions/reference/excursions.inla.md),
except that the `iNIQC` method is not available. With the `NI` and
`NIQC` methods, the joint probability of a region is the mixture over
the hyperparameter configurations of INLA, and the regions are grown
using the correlations of the configuration with the largest posterior
density.

Models fitted with `inlabru` are handled in the same way as in
[`excursions.inla()`](https://davidbolin.github.io/excursions/reference/excursions.inla.md).
To compute regions for the linear predictor at a set of locations, such
as the nodes of a mesh, add a likelihood component with `NA`
observations at the locations and a tag, and use `name = "APredictor"`
and `ind = bru_index(result, tag)`. The graph is then given for the
locations, for example as the mesh.

## Note

This function requires the `INLA` package, which is not a CRAN package.
See <https://www.r-inla.org/download-install> for easy installation
instructions.

## See also

[`excursions.regions()`](https://davidbolin.github.io/excursions/reference/excursions.regions.md),
[`excursions.inla()`](https://davidbolin.github.io/excursions/reference/excursions.inla.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
if (require.nowarnings("INLA") && require.nowarnings("inlabru")) {
  ## Simulate data on a mesh
  x <- seq(from = 0, to = 10, length.out = 20)
  lattice <- fmesher::fm_lattice_2d(x = x, y = x)
  mesh <- fmesher::fm_rcdt_2d_inla(
    lattice = lattice, extend = FALSE, refine = FALSE
  )
  Q <- fmesher::fm_matern_precision(mesh, alpha = 2, rho = 3, sigma = 1)
  field <- fmesher::fm_sample(n = 1, Q = Q)
  obs.loc <- matrix(runif(200) * 10, 100, 2)
  y <- as.vector(fmesher::fm_basis(mesh, loc = obs.loc) %*% field) +
    rnorm(100) * 0.3

  ## Fit the model with inlabru, with NA observations at the mesh nodes
  matern <- INLA::inla.spde2.pcmatern(mesh,
    prior.range = c(1, 0.5), prior.sigma = c(1, 0.5)
  )
  data <- data.frame(x1 = obs.loc[, 1], x2 = obs.loc[, 2], y = y)
  data.prd <- data.frame(x1 = mesh$loc[, 1], x2 = mesh$loc[, 2], y = NA)
  fit <- inlabru::bru(
    ~ Intercept(1) + field(cbind(x1, x2), model = matern),
    inlabru::bru_obs(y ~ ., family = "normal", data = data),
    inlabru::bru_obs(y ~ ., family = "normal", data = data.prd, tag = "prd"),
    options = list(control.compute = list(return.marginals.predictor = TRUE))
  )

  ## Connected regions where the field exceeds 0
  res <- excursions.regions.inla(fit,
    name = "APredictor", ind = inlabru::bru_index(fit, "prd"),
    graph = mesh, alpha = 0.1, u = 0, type = ">", method = "QC",
    min.size = 5, prune.ind = TRUE
  )
  lengths(res$regions)
}
} # }
```
