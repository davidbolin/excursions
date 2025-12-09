# Contour maps and contour map quality measures using Monte Carlo samples

`contourmap.mc` is used for calculating contour maps and quality
measures for contour maps based on Monte Carlo samples of a model.

## Usage

``` r
contourmap.mc(
  samples,
  n.levels,
  ind,
  levels,
  type = c("standard", "equalarea", "P0-optimal", "P1-optimal", "P2-optimal"),
  compute = list(F = TRUE, measures = NULL),
  alpha,
  verbose = FALSE
)
```

## Arguments

- samples:

  Matrix with model Monte Carlo samples. Each column contains a sample
  of the model.

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

- alpha:

  Maximal error probability in contour map function (default=0.1).

- verbose:

  Set to TRUE for verbose mode (optional).

## Value

`contourmap` returns an object of class "excurobj" with the following
elements

- u :

  Contour levels used in the contour map.

- n.levels :

  The number of contours used.

- u.e :

  The values associated with the level sets `G_k`.

- G :

  A vector which shows which of the level sets `G_k` each node belongs
  to.

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

The contour map is computed for the empirical mean of the samples. See
[`contourmap()`](https://davidbolin.github.io/excursions/reference/contourmap.md)
and
[`contourmap.inla()`](https://davidbolin.github.io/excursions/reference/contourmap.inla.md)
for further details.

## References

Bolin, D. and Lindgren, F. (2017) *Quantifying the uncertainty of
contour maps*, Journal of Computational and Graphical Statistics, 26:3,
513-524.

Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion
Sets and Related Quantities Using excursions*, Journal of Statistical
Software, 86(5), 1–20.

## See also

[`contourmap()`](https://davidbolin.github.io/excursions/reference/contourmap.md),
[`contourmap.inla()`](https://davidbolin.github.io/excursions/reference/contourmap.inla.md),
[`contourmap.colors()`](https://davidbolin.github.io/excursions/reference/contourmap.colors.md)

## Author

David Bolin <davidbolin@gmail.com>

## Examples

``` r
n <- 100
Q <- Matrix(toeplitz(c(1, -0.5, rep(0, n - 2))))
mu <- seq(-5, 5, length = n)
## Sample the model 100 times (increase for better estimate)
X <- mu + solve(chol(Q), matrix(rnorm(n = n * 100), nrow = n, ncol = 100))

lp <- contourmap.mc(X, n.levels = 2, compute = list(F = FALSE, measures = c("P1", "P2")))

# plot contourmap
plot(lp$map)

# display quality measures
c(lp$P1, lp$P2)
#> [1] 0.2 0.0
```
