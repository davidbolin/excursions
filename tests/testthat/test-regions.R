regions.testdata <- function() {
  n <- 40
  phi <- 0.9
  Q <- sparseMatrix(
    i = c(1:n, 2:n), j = c(1:n, 1:(n - 1)),
    x = c(1, rep(1 + phi^2, n - 2), 1, rep(-phi, n - 1)) / (1 - phi^2),
    dims = c(n, n), symmetric = TRUE
  )
  x <- seq(0, 1, length.out = n)
  mu <- 4 * exp(-(x - 0.3)^2 / 0.01) + 3 * exp(-(x - 0.75)^2 / 0.005) - 1
  list(n = n, Q = Q, mu = mu)
}

test_that("Bivariate normal probabilities", {
  ref <- function(h, k, r) {
    integrate(function(x) dnorm(x) * pnorm((k - r * x) / sqrt(1 - r^2)),
      -Inf, h,
      rel.tol = 1e-12, abs.tol = 1e-15
    )$value
  }
  g <- expand.grid(
    h = c(-3, -0.4, 0.7, 2.5), k = c(-2.2, 0.9, 3),
    r = c(-0.99, -0.95, -0.5, 0, 0.6, 0.93, 0.98)
  )
  p <- .Call("regions_bvn_lower", g$h, g$k, g$r, PACKAGE = "excursions")
  expect_equal(p, mapply(ref, g$h, g$k, g$r), tolerance = 1e-10)
})

test_that("Selected inverse", {
  m <- 6
  D <- bandSparse(m, k = -1:1, diagonals = list(
    rep(-1, m - 1), rep(2.1, m), rep(-1, m - 1)
  ))
  Q <- kronecker(Diagonal(m), D) + kronecker(D, Diagonal(m))
  Q <- Q %*% Q
  S <- solve(as.matrix(Q))
  sel <- private.selected.inverse(Q)
  expect_equal(sel$vars, diag(S), tolerance = 1e-10)
  tr <- summary(as(Q, "generalMatrix"))
  expect_equal(sel$cov(tr$i, tr$j), S[cbind(tr$i, tr$j)], tolerance = 1e-10)
})

test_that("Prominence of local maxima", {
  ## A path with maxima at nodes 2, 4 and 6. Node 6 is the highest, node 2
  ## reaches it through the minimum 0.5 at node 5, and node 4 reaches node 2
  ## through the minimum 1 at node 3.
  z <- c(0, 3, 1, 2, 0.5, 4, 0)
  n <- length(z)
  G <- private.regions.graph(bandSparse(n, k = 1, symmetric = TRUE), NULL, n)
  prom <- .Call("regions_prominence", G@p, G@i, z, as.integer(order(-z)),
    PACKAGE = "excursions"
  )
  expect_equal(prom, c(0, 2.5, 0, 1, 0, Inf, 0))
  ## Within a subset of the nodes, the components are separate
  prom <- .Call("regions_prominence", G@p, G@i, z, as.integer(c(2, 4, 3, 1)),
    PACKAGE = "excursions"
  )
  expect_equal(prom, c(0, Inf, 0, 1, 0, 0, 0))
})

test_that("Regions are connected, disjoint and satisfy the probability", {
  data <- regions.testdata()
  alpha <- 0.1
  res <- excursions.regions(
    alpha = alpha, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = 1, max.threads = 1
  )
  expect_true(length(res$regions) >= 2)
  for (R in res$regions) {
    expect_true(all(diff(R) == 1))
  }
  expect_false(anyDuplicated(unlist(res$regions)) > 0)
  expect_true(all(res$P >= 1 - alpha))
  expect_true(all(res$rho[unlist(res$regions)] >= 1 - alpha))
  expect_equal(lengths(res$regions), sort(lengths(res$regions), decreasing = TRUE))
  expect_equal(which(res$labels == 1), res$regions[[1]])
  expect_equal(which(res$E == 1), res$regions[[1]])

  ## The probability of the largest region agrees with gaussint
  R <- res$regions[[1]]
  a <- rep(-Inf, data$n)
  a[R] <- 0
  p <- gaussint(
    mu = data$mu, Q = data$Q, a = a, b = rep(Inf, data$n),
    use.reordering = "limits", seed = 2, max.threads = 1
  )$P
  expect_equal(res$P[1], p, tolerance = 0.01)

  ## The largest region is at least as large as the largest connected part
  ## of the excursion set
  ex <- excursions(
    alpha = alpha, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = 1, max.threads = 1
  )
  runs <- rle(ex$E == 1)
  expect_true(length(R) >= max(runs$lengths[runs$values]))
})

test_that("Regions for type < mirror type >", {
  data <- regions.testdata()
  r1 <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = 1, max.threads = 1
  )
  r2 <- excursions.regions(
    alpha = 0.1, u = 0, mu = -data$mu, Q = data$Q, type = "<",
    seed = 1, max.threads = 1
  )
  expect_equal(r1$regions, r2$regions)
  ## The sampler uses the random numbers differently for lower and upper
  ## limits, so the probabilities only agree up to Monte Carlo error
  expect_equal(r1$P, r2$P, tolerance = 0.01)
})

test_that("Region options", {
  data <- regions.testdata()
  res <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    min.size = 3, max.regions = 1, seed = 1, max.threads = 1
  )
  expect_length(res$regions, 1)

  res <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    min.size = 3, seed = 1, max.threads = 1
  )
  expect_true(all(lengths(res$regions) >= 3))

  ## Without edges, each candidate node is its own region
  res <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    graph = Diagonal(data$n), seed = 1, max.threads = 1
  )
  expect_true(all(lengths(res$regions) == 1))
  expect_equal(sort(unlist(res$regions)), which(res$rho >= 0.9))

  ## Growth by marginal probabilities, restricted to a set of nodes
  res <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    growth = "rho", ind = 1:20, seed = 1, max.threads = 1
  )
  expect_true(all(unlist(res$regions) <= 20))

  ## The first start point is the same, so more start points cannot give a
  ## smaller region
  r1 <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    n.starts = 1, max.regions = 1, seed = 1, max.threads = 1
  )
  r5 <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    n.starts = 5, max.regions = 1, seed = 1, max.threads = 1
  )
  expect_true(length(r5$regions[[1]]) >= length(r1$regions[[1]]))

  ## With a large minimum prominence, only the largest maximum is used
  r.prom <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    n.starts = 5, min.prominence = 100, seed = 1, max.threads = 1
  )
  r.one <- excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    n.starts = 1, seed = 1, max.threads = 1
  )
  expect_equal(r.prom$regions, r.one$regions)

  expect_error(excursions.regions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = "="
  ))
})

test_that("Excursion functions of the regions", {
  data <- regions.testdata()
  alpha <- 0.1
  res <- excursions.regions(
    alpha = alpha, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = 1, max.threads = 1
  )
  expect_s3_class(res, "excurobj")
  expect_output(print(res), "Connected positive excursion regions")
  expect_equal(dim(res$F), c(data$n, length(res$regions)))
  for (r in seq_along(res$regions)) {
    R <- res$regions[[r]]
    Fr <- as.vector(res$F[, r])
    ## The region is where the excursion function is at least 1 - alpha,
    ## and the excursion function decreases along the growth order
    expect_equal(which(Fr >= 1 - alpha), R)
    expect_equal(min(Fr[R]), res$P[r])
    ## It is zero at the nodes of the other regions, and given for the
    ## neighbours of the region
    expect_true(all(Fr[res$labels > 0 & res$labels != r] == 0))
    nb <- setdiff(c(min(R) - 1, max(R) + 1), c(0, data$n + 1))
    nb <- nb[res$labels[nb] == 0]
    expect_true(all(Fr[nb] > 0 & Fr[nb] < 1 - alpha))
  }

  ## The value at a neighbour is the joint probability of the region and the
  ## neighbour
  R <- res$regions[[1]]
  j <- max(R) + 1
  a <- rep(-Inf, data$n)
  a[c(R, j)] <- 0
  p <- gaussint(
    mu = data$mu, Q = data$Q, a = a, b = rep(Inf, data$n),
    use.reordering = "limits", seed = 2, max.threads = 1
  )$P
  expect_equal(res$F[j, 1], p, tolerance = 0.01)
})

test_that("Continuous domain regions", {
  skip_if_not_installed("sp")
  x <- seq(0, 10, length.out = 12)
  mesh <- fmesher::fm_rcdt_2d_inla(
    lattice = fmesher::fm_lattice_2d(x = x, y = x),
    extend = FALSE, refine = FALSE
  )
  Q <- fmesher::fm_matern_precision(mesh, alpha = 2, rho = 2, sigma = 1)
  mu <- 2.5 * cos(mesh$loc[, 1] / 1.5) * cos(mesh$loc[, 2] / 2)
  res <- excursions.regions(
    alpha = 0.1, u = 0, mu = mu, Q = Q, type = ">", graph = mesh,
    min.size = 2, seed = 1, max.threads = 1
  )
  expect_true(length(res$regions) >= 2)
  sets <- continuous(res, mesh)
  expect_equal(
    vapply(sets$M@polygons, function(p) p@ID, ""),
    as.character(seq_along(res$regions))
  )
  expect_equal(ncol(sets$F), length(res$regions))
  ## Each region contains its nodes, and no nodes of the other regions
  loc <- sp::SpatialPoints(mesh$loc[, 1:2])
  for (r in seq_along(res$regions)) {
    inside <- !is.na(sp::over(loc, sets$M[as.character(r)]))
    expect_true(all(inside[res$labels == r]))
    expect_false(any(inside[res$labels > 0 & res$labels != r]))
  }
  sets.fm <- continuous(res, mesh, output = "fm")
  expect_equal(sort(unique(sets.fm$M$grp)), seq_along(res$regions))
})
