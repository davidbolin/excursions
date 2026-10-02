test_that("Contourmap, test ind", {
  data <- integration.testdata1()
  ind1 <- c(1, 2, 3, 4)
  ind2 <- c(4, 3, 2, 1)
  ind3 <- rep(FALSE, data$n)
  ind3[1:4] <- TRUE

  res1 <- contourmap(data$mu, data$Q,
    n.levels = 2, ind = ind1,
    seed = data$seed, alpha = 0.1, max.threads = 1
  )
  res2 <- contourmap(data$mu, data$Q,
    n.levels = 2, ind = ind2,
    seed = data$seed, alpha = 0.1, max.threads = 1
  )
  res3 <- contourmap(data$mu, data$Q,
    n.levels = 2, ind = ind3,
    seed = data$seed, alpha = 0.1, max.threads = 1
  )

  expect_equal(res1$F, res2$F, tolerance = 1e-7)
  expect_equal(res2$F, res3$F, tolerance = 1e-7)
})


test_that("Contourmap, P measures", {
  data <- integration.testdata1()

  res1 <- contourmap(data$mu, data$Q,
    n.levels = 4,
    seed = data$seed, alpha = 0.1, max.threads = 1,
    compute = list(F = FALSE, measures = c("P2", "P1"))
  )

  expect_equal(res1$P1, 0.9217417, tolerance = 1e-3)
  expect_equal(res1$P2, 0.405841, tolerance = 1e-3)
})


test_that("Contour maps are unchanged", {
  d <- testdata.spde(10)
  r <- contourmap(d$mu, d$Q,
    n.levels = 3, seed = d$seed, alpha = 0.1, max.threads = 1,
    n.iter = 2000, size.tol = NULL,
    compute = list(F = TRUE, measures = c(
      "P0", "P1", "P2",
      "P0-bound", "P1-bound", "P2-bound"
    ))
  )
  expect_equal(
    c(
      F = ref.summary(r$F), E = ref.summary(r$E), M = ref.summary(r$M),
      P0 = r$P0, P1 = r$P1, P2 = r$P2,
      P0b = r$P0.bound, P1b = r$P1.bound, P2b = r$P2.bound
    ),
    REF$cm,
    tolerance = 1e-8
  )
})

test_that("Contour maps compute the variances once", {
  d <- testdata.spde(8)
  n.calls <- 0
  variances <- excursions.variances
  local_mocked_bindings(excursions.variances = function(...) {
    n.calls <<- n.calls + 1
    variances(...)
  })
  r <- contourmap(d$mu, d$Q,
    n.levels = 3, seed = d$seed, alpha = 0.1, max.threads = 1,
    n.iter = 500,
    compute = list(F = TRUE, measures = c("P0", "P0-bound", "P1-bound"))
  )
  expect_equal(n.calls, 1)

  ## The results are the same as when the variances are supplied
  vars <- variances(Q = d$Q)
  r.vars <- contourmap(d$mu, d$Q,
    vars = vars,
    n.levels = 3, seed = d$seed, alpha = 0.1, max.threads = 1,
    n.iter = 500,
    compute = list(F = TRUE, measures = c("P0", "P0-bound", "P1-bound"))
  )
  expect_equal(r$F, r.vars$F)
  expect_equal(r$P0.bound, r.vars$P0.bound)
  expect_equal(r$P1.bound, r.vars$P1.bound)
})

test_that("Contour map function with a Cholesky factor matches the precision matrix", {
  d <- testdata.spde(8)
  lp <- contourmap(d$mu, d$Q, n.levels = 2, seed = d$seed, alpha = 0.1, max.threads = 1, n.iter = 500)
  f.Q <- excursions:::contourfunction(
    lp = lp, mu = d$mu, Q = d$Q, alpha = 0.1, F.limit = 0.1,
    seed = d$seed, max.threads = 1, n.iter = 500
  )
  f.L <- excursions:::contourfunction(
    lp = lp, mu = d$mu, Q.chol = chol(d$Q), alpha = 0.1, F.limit = 0.1,
    seed = d$seed, max.threads = 1, n.iter = 500
  )
  expect_equal(f.L$F, f.Q$F, tolerance = 1e-10)
})

test_that("contourmap with an adaptive number of iterations", {
  d <- testdata.spde(10)
  run <- function(...) {
    contourmap(d$mu, d$Q,
      n.levels = 2, alpha = 0.1,
      compute = list(F = TRUE, measures = c("P1", "P2")),
      seed = d$seed, max.threads = 1, ...
    )
  }
  ## With a large tol, only the first batch of 1000 iterations is used
  r0 <- run(n.iter = 1000)
  r1 <- run(n.iter = 10000, tol = 1)
  expect_identical(r1$F, r0$F)
  expect_identical(r1$P1, r0$P1)
  expect_identical(r1$P2, r0$P2)
  expect_equal(r1$meta$tol, 1)
})
