test_that("Continous on contourmap, R2 mesh", {
  skip_on_cran()
  skip_if_not_installed("sp")

  data <- integration.testdata1()
  res1 <- contourmap(data$mu, data$Q,
    n.levels = 2,
    seed = data$seed, alpha = 0.1, max.threads = 1,
    compute = list(
      F = TRUE,
      measures = c("P2", "P1", "P0")
    )
  )

  bnd <- fmesher::fm_segm(
    cbind(
      c(0, 1, 2, 2, 2, 1, 0, 0),
      c(0, 0, 0, 1, 2, 2, 2, 1)
    ),
    is.bnd = TRUE
  )
  loc <- cbind(
    c(0.5, 1.5, 1),
    c(0.5, 0.5, 1.5)
  )
  mesh <- fmesher::fm_rcdt_2d_inla(loc = loc, boundary = bnd)
  expect_equal(data$n, mesh$n)

  GQ_direct <- gaussquad(mesh, method = "direct")
  GQ_make_A <- gaussquad(mesh, method = "make.A")
  expect_equal(GQ_direct$A, GQ_make_A$A, tolerance = 1e-15)

  res2 <- continuous(res1, mesh, method = "linear")

  expect_error(
    continuous(res1, mesh, method = "linear", output = "fm"),
    "Output format 'fm' not supported for 'calc.credible = TRUE'."
  )
  res3 <- continuous(res1, mesh,
    method = "linear", output = "fm",
    calc.credible = FALSE
  )
  res4 <- continuous(res1, mesh,
    method = "log", output = "fm",
    calc.credible = FALSE
  )
  expect_s4_class(res2$M, "SpatialPolygons")
  expect_s3_class(res3$M, "fm_segm")
  expect_s3_class(res4$M, "fm_segm")

  testthat::skip_if_not_installed("fmesher", "0.1.2.9003")
  # 0.1.2 had a bug for empty segments.
  res5 <- continuous(res1, mesh,
    method = "step", output = "fm",
    calc.credible = FALSE
  )
  expect_s3_class(res5$M, "fm_segm")
})

test_that("Continous on contourmap, M mesh", {
  skip_on_cran()
  skip_if_not_installed("sp")

  data <- integration.testdata1()
  res1 <- contourmap(data$mu, data$Q,
    n.levels = 2,
    seed = data$seed, alpha = 0.1, max.threads = 1
  )

  bnd <- fmesher::fm_segm(
    cbind(
      c(0, 1, 2, 2, 2, 1, 0, 0),
      c(0, 0, 0, 1, 2, 2, 2, 1)
    ),
    is.bnd = TRUE
  )
  loc <- cbind(
    c(0.5, 1.5, 1),
    c(0.5, 0.5, 1.5)
  )
  mesh <- fmesher::fm_rcdt_2d_inla(loc = loc, boundary = bnd)
  expect_equal(data$n, mesh$n)

  # Alter z-coordinates and mark as general manifold
  mesh$loc[, 3] <- seq_len(mesh$n)
  mesh$manifold <- "M2"

  res2 <- continuous(res1, mesh,
    method = "linear",
    output = "fm",
    calc.credible = FALSE
  )

  expect_s3_class(res2$M, "fm_segm")

  res3 <- continuous(res1, mesh,
    method = "log",
    output = "fm",
    calc.credible = FALSE
  )

  expect_s3_class(res3$M, "fm_segm")

  testthat::skip_if_not_installed("fmesher", "0.1.2.9003")
  res4 <- continuous(res1, mesh,
    method = "step",
    output = "fm",
    calc.credible = FALSE
  )
  expect_s3_class(res4$M, "fm_segm")
})

test_that("Continuous on excursions, lattice with partial ind", {
  skip_on_cran()
  skip_if_not_installed("sp")

  nxy <- 10
  x <- seq(0, 1, length.out = nxy)
  lattice <- fmesher::fm_lattice_2d(x = x, y = x)
  mesh <- fmesher::fm_rcdt_2d_inla(lattice = lattice, extend = FALSE, refine = FALSE)
  Q <- fmesher::fm_matern_precision(mesh, alpha = 2, rho = 0.3, sigma = 1)
  # Order mean and precision as the lattice nodes
  reo <- mesh$idx$lattice
  Q <- Q[reo, reo]
  mu <- 3 * (lattice$loc[, 1] + lattice$loc[, 2] - 1)

  # Leave out a corner block, so that the active set is a strict subset
  ind <- which(!(lattice$loc[, 1] < 0.35 & lattice$loc[, 2] < 0.35))
  ex <- excursions(
    alpha = 0.1, u = 0, mu = mu, Q = Q, type = ">",
    ind = ind, F.limit = 1, seed = 1:6, max.threads = 1
  )
  expect_true(any(ex$E[ind] == 1))

  res <- continuous(ex, lattice, alpha = 0.1, method = "linear")
  expect_s4_class(res$M, "SpatialPolygons")

  # The interpolated F must equal the input F at the active lattice nodes
  loc.key <- function(loc) paste(round(loc[, 1], 8), round(loc[, 2], 8))
  vtx <- match(loc.key(lattice$loc[ind, , drop = FALSE]), loc.key(res$F.geometry$loc))
  expect_false(anyNA(vtx))
  F.in <- ex$F[ind]
  F.in[is.na(F.in)] <- 0
  expect_equal(res$F[vtx], F.in, tolerance = 1e-10)
})
