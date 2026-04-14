# Contour maps and contour map quality measures for latent Gaussian models

`contourmap` is used for calculating contour maps and quality measures
for contour maps for Gaussian models.

## Usage

``` r
contourmap(
  mu,
  Q,
  vars,
  n.levels,
  ind,
  levels,
  type = c("standard", "pretty", "equalarea", "P0-optimal", "P1-optimal", "P2-optimal"),
  compute = list(F = TRUE, measures = NULL),
  use.marginals = TRUE,
  alpha,
  F.limit,
  n.iter = 10000,
  verbose = FALSE,
  max.threads = 0,
  seed = NULL
)
```

## Arguments

- mu:

  Expectation vector.

- Q:

  Precision matrix.

- vars:

  Precomputed marginal variances (optional).

- n.levels:

  Number of levels in contour map.

- ind:

  Indices of the nodes that should be analyzed (optional).

- levels:

  Levels to use in contour map.

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

  'P0-optimal'

  :   Levels chosen to maximize the P0 measure.

  'P1-optimal'

  :   Levels chosen to maximize the P1 measure.

  'P2-optimal'

  :   Levels chosen to maximize the P2 measure.

- compute:

  A list with quality indices to compute

  'F':

  :   TRUE/FALSE indicating whether the contour map function should be
      computed (default TRUE).

  'measures':

  :   A list with the quality measures to compute ("P0", "P1", "P2") or
      corresponding bounds based only on the marginal probabilities
      ("P0-bound", "P1-bound", "P2-bound").

- use.marginals:

  Only marginal distributions are used when finding P-optimal maps
  (default TRUE).

- alpha:

  Maximal error probability in contour map function (default=1).

- F.limit:

  The limit value for the computation of the F function. F is set to NA
  for all nodes where `F < 1-F.limit`. Default is `F.limit = alpha`.

- n.iter:

  Number or iterations in the MC sampler that is used for calculating
  the quantities in `compute`. The default value is 10000.

- verbose:

  Set to TRUE for verbose mode (optional).

- max.threads:

  Decides the number of threads the program can use. Set to 0 for using
  the maximum number of threads allowed by the system (default).

- seed:

  Random seed (optional).

## Value

`contourmap` returns an object of class "excurobj" with the following
elements

- u :

  Contour levels used in the contour map.

- n.levels :

  The number of contours used.

- u.e :

  The values associated with the level sets G_k.

- G :

  A vector which shows which of the level sets G_k each node belongs to.

- map :

  Representation of the contour map with `map[i]=u.e[k]` if i is in
  `G_k`.

- F :

  The contour map function (if computed).

- M :

  Contour avoiding sets (if `F` is computed). \\M=-1\\ for all
  non-significant nodes and \\M=k\\ for nodes that belong to \\M_k\\.

- P0/P1/P2 :

  Calculated quality measures (if computed).

- P0bound/P1bound/P2bound :

  Calculated upper bounds quality measures (if computed).

- meta :

  A list containing various information about the calculation.

## Details

The Gaussian model is specified using the mean `mu` and the precision
matrix `Q`. The contour map is then computed for the mean, using either
the contour levels specified in `levels`, or `n.levels` contours that
are placed according to the argument `type`.

A number of quality measures can be computed based based on the
specified contour map and the Gaussian distribution. What should be
computed is specified using the `compute` argument. For details on these
quanties, see the reference below.

## References

Bolin, D. and Lindgren, F. (2017) *Quantifying the uncertainty of
contour maps*, Journal of Computational and Graphical Statistics, vol
26, no 3, pp 513-524.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, vol 86, no 1, pp 1-20.

## See also

[`contourmap.inla()`](https://davidbolin.github.io/excursions/reference/contourmap.inla.md),
[`contourmap.mc()`](https://davidbolin.github.io/excursions/reference/contourmap.mc.md),
[`contourmap.colors()`](https://davidbolin.github.io/excursions/reference/contourmap.colors.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
n <- 10
Q <- Matrix(toeplitz(c(1, -0.5, rep(0, n - 2))))
mu <- seq(-5, 5, length = n)
lp <- contourmap(mu, Q,
  n.levels = 2,
  compute = list(F = FALSE, measures = c("P1", "P2")),
  max.threads = 1
)
# Plot the contourmap
plot(lp$map)

# Display the quality measures
cat(c(lp$P1, lp$P2))
#> 0.908944 0.4677515
```
