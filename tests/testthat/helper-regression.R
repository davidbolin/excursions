## Test problems shared by the regression tests. The reference values in the
## tests were computed with these problems, so do not change them without
## recomputing the references.

## Precision matrix of a second order SPDE on an m x m lattice
testdata.spde <- function(m = 15, kappa = 0.5) {
  D <- Matrix::bandSparse(m,
    k = c(-1, 0, 1),
    diagonals = list(rep(-1, m - 1), rep(2, m), rep(-1, m - 1))
  )
  I <- Matrix::Diagonal(m)
  K <- kappa^2 * kronecker(I, I) + kronecker(D, I) + kronecker(I, D)
  Q <- Matrix::forceSymmetric(Matrix::crossprod(K))
  x <- seq(0, 1, length.out = m)
  mu <- 2 * as.vector(outer(sin(x * 6), cos(x * 5)))
  list(Q = Q, mu = mu, n = m^2, m = m, x = x, seed = 1:6)
}

## Samples from the SPDE model, one column per sample
testdata.spde.samples <- function(data, n.samples = 500) {
  L <- Matrix::chol(data$Q)
  Z <- withr::with_seed(1, matrix(stats::rnorm(data$n * n.samples), data$n))
  data$mu + as.matrix(Matrix::solve(L, Z))
}

## Gaussian mixture of three scaled versions of the SPDE model
testdata.mixture <- function(m = 8) {
  data <- testdata.spde(m)
  list(
    mu = list(data$mu, data$mu + 0.2, data$mu - 0.1),
    Q = list(data$Q, 1.3 * data$Q, 0.7 * data$Q),
    w = c(0.5, 0.3, 0.2),
    n = data$n,
    seed = 1:6
  )
}

## Triangulation of the unit square with a field that has several contours
testdata.mesh <- function(m = 25) {
  mesh <- fmesher::fm_rcdt_2d_inla(
    globe = NULL,
    lattice = fmesher::fm_lattice_2d(
      x = seq(0, 1, length.out = m),
      y = seq(0, 1, length.out = m)
    )
  )
  z <- sin(mesh$loc[, 1] * 9) + cos(mesh$loc[, 2] * 7) + 0.3 * mesh$loc[, 1]
  list(mesh = mesh, z = z)
}

## Weighted column sums of the rows of a matrix after sorting them, which
## does not depend on the order of the rows. The values are rounded before
## sorting, so that rounding differences between platforms cannot change the
## order of rows with (almost) equal values.
sorted.summary <- function(M, digits = 8) {
  M <- as.matrix(M)
  if (nrow(M) == 0) {
    return(rep(0, 2 * ncol(M)))
  }
  R <- round(M, digits)
  o <- do.call(order, lapply(seq_len(ncol(R)), function(k) R[, k]))
  c(colSums(M), colSums(M[o, , drop = FALSE] * seq_along(o)))
}

## Summary of the output of tricontour, for comparisons with reference values
## without storing the whole output. The numbering of the vertices and the
## order of the edges depend on the numbering of the mesh, which can differ
## between platforms, so the summary only uses the geometry: the vertex
## coordinates, and the edges as the coordinates of their end points together
## with their group.
fingerprint.tricontour <- function(x) {
  loc <- x$loc[, 1:2, drop = FALSE]
  edges <- cbind(loc[x$idx[, 1], , drop = FALSE], loc[x$idx[, 2], , drop = FALSE], x$grp)
  list(
    nloc = nrow(loc),
    loc = sorted.summary(loc),
    nidx = nrow(x$idx),
    edges = sorted.summary(edges),
    grp = tabulate(x$grp)
  )
}

## Summary of the output of connect.segments, for segments given by node
## indices in a fixed order
fingerprint.segments <- function(x) {
  list(
    n = length(x$sequences),
    len = vapply(x$sequences, length, 1L),
    seq = vapply(x$sequences, function(s) sum(s * seq_along(s)), 0),
    seg = vapply(x$seg, function(s) sum(s * seq_along(s)), 0),
    grp = vapply(x$grp, function(s) sum(s * seq_along(s)), 0)
  )
}

## Summary of the output of connect.segments that only depends on the
## geometry: for each sequence, its length, whether it is closed, its groups,
## and the coordinates of its nodes (without the repeated first node of a
## closed sequence, whose choice depends on the numbering).
fingerprint.segments.geometry <- function(x, loc) {
  seqs <- t(vapply(seq_along(x$sequences), function(k) {
    s <- x$sequences[[k]]
    closed <- s[1] == s[length(s)]
    nodes <- if (closed) s[-length(s)] else s
    c(
      length(s), closed, sum(x$grp[[k]]),
      colSums(loc[nodes, 1:2, drop = FALSE])
    )
  }, numeric(5)))
  list(n = length(x$sequences), seqs = sorted.summary(seqs))
}

## Summary of the output of continuous(): the excursion function at the
## nodes of the geometry, paired with their coordinates, and the rings of the
## polygons, without their start point.
fingerprint.continuous <- function(r) {
  rings <- do.call(rbind, lapply(r$M@polygons, function(p) {
    t(vapply(p@Polygons, function(q) {
      crd <- q@coords[-nrow(q@coords), , drop = FALSE]
      c(q@hole, nrow(crd), q@area, colSums(crd))
    }, numeric(5)))
  }))
  F <- r$F
  F[is.na(F)] <- -1
  list(
    F = sorted.summary(cbind(r$F.geometry$loc[, 1:2], F)),
    nrings = nrow(rings),
    rings = sorted.summary(rings)
  )
}

## Sum and weighted sum of a vector, with missing values set to -1
ref.summary <- function(v) {
  v[is.na(v)] <- -1
  c(sum = sum(v), wsum = sum(v * seq_along(v)))
}

## Check a simultaneous confidence band for the mixture in testdata.mixture()
check.mixture.band <- function(r, d, alpha = 0.1) {
  n <- d$n
  testthat::expect_length(r$a, n)
  testthat::expect_true(all(is.finite(r$a)) && all(is.finite(r$b)))
  ## The simultaneous band contains the marginal band
  testthat::expect_true(all(r$a <= r$a.marginal & r$b >= r$b.marginal))
  ## Mean and variance of the mixture
  vars <- lapply(d$Q, function(Q) diag(solve(Q)))
  mix.mean <- Reduce(`+`, Map(`*`, d$w, d$mu))
  mix.vars <- Reduce(`+`, Map(function(w, m, v) w * (v + m^2), d$w, d$mu, vars)) -
    mix.mean^2
  testthat::expect_equal(r$mean, mix.mean, tolerance = 1e-10)
  testthat::expect_equal(r$vars, mix.vars, tolerance = 1e-10)
  ## Empirical coverage for new samples from the mixture
  X <- withr::with_seed(42, {
    k <- sample(seq_along(d$w), 4000, replace = TRUE, prob = d$w)
    vapply(k, function(kk) {
      d$mu[[kk]] + as.vector(Matrix::solve(
        Matrix::chol(d$Q[[kk]]),
        stats::rnorm(n)
      ))
    }, numeric(n))
  })
  coverage <- mean(colSums(X < r$a | X > r$b) == 0)
  testthat::expect_equal(coverage, 1 - alpha, tolerance = 0.03 / (1 - alpha))
}
