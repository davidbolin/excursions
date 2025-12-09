# Extract a part of a mesh

Extracts a part of a mesh

## Usage

``` r
submesh.mesh(z, mesh)
```

## Arguments

- z:

  A matrix with values indicating which nodes that should be present in
  the submesh.

- mesh:

  An `fm_mesh_2d` object.

## Value

An `fm_mesh_2d` object.

## Author

Finn Lindgren <finn.lindgren@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
if (require(fmesher)) {
  nxy <- 30
  x <- seq(from = 0, to = 4, length.out = nxy)
  lattice <- fm_lattice_2d(x = x, y = x)
  mesh <- fm_rcdt_2d_inla(lattice = lattice, extend = FALSE, refine = FALSE)

  # extract a part of the mesh inside a circle
  xy.in <- rowSums((mesh$loc[, 1:2] - 2)^2) < 1
  submesh <- submesh.mesh(matrix(xy.in, nxy, nxy), mesh)
  plot(mesh$loc[, 1:2], pch = 20)
  lines(2 + cos(seq(0, 2 * pi, length.out = 100)), 2 + sin(seq(0, 2 * pi, length.out = 100)))
  plot(submesh, add = TRUE)
  points(mesh$loc[xy.in, 1:2], col = "2")
}
} # }
```
