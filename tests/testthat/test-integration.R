test_that("Integration L", {
  data <- integration.testdata1()
  prob1 <- gaussint(
    Q.chol = data$L, a = data$a, b = data$b,
    seed = data$seed, max.threads = 1
  )
  expect_equal(prob1$P[1], 0.9680023, tolerance = 1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance = 1e-6)
})

test_that("Integration Q", {
  data <- integration.testdata1()
  prob1 <- gaussint(
    Q = data$Q, a = data$a, b = data$b,
    seed = data$seed, max.threads = 1
  )
  expect_equal(prob1$P[1], 0.9680023, tolerance = 1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance = 1e-6)
})

test_that("Integration mu", {
  data <- integration.testdata1()
  prob1 <- gaussint(
    Q = data$Q, mu = data$mu, a = data$a + data$mu,
    b = data$b + data$mu, seed = data$seed,
    max.threads = 1
  )
  expect_equal(prob1$P[1], 0.9680023, tolerance = 1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance = 1e-6)
})

test_that("Integration limit", {
  data <- integration.testdata1()
  prob1 <- gaussint(
    Q = data$Q, a = data$a, b = data$b,
    seed = data$seed, lim = 0.97,
    max.threads = 1
  )

  prob2 <- gaussint(
    Q = data$Q, a = data$a, b = data$b,
    seed = data$seed, lim = 0.9,
    max.threads = 1
  )

  expect_equal(prob1$P[1], 0.0, tolerance = 1e-7)
  expect_equal(prob2$P[1], 0.9680023, tolerance = 1e-6)
})

test_that("Integration reordering", {
  data <- integration.testdata1()
  prob1 <- gaussint(
    Q = data$Q, mu = data$mu, a = data$a + data$mu,
    b = data$b + data$mu, seed = data$seed,
    max.threads = 1, use.reordering = "sparsity"
  )
  expect_equal(prob1$P[1], 0.9680023, tolerance = 1e-5)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance = 1e-5)
})


## Summary of a gaussint result, as stored in REF
gaussint.ref <- function(...) {
  r <- gaussint(..., n.iter = 2000, max.threads = 1)
  c(P = r$P, E = r$E, Pv = ref.summary(r$Pv), Ev = ref.summary(r$Ev))
}

test_that("Integration results are unchanged", {
  d <- testdata.spde(10)
  a <- d$mu - 2.5
  b <- d$mu + 2.5
  lind <- seq_len(d$n) %in% 30:70
  tol <- 1e-8

  expect_equal(gaussint.ref(Q = d$Q, a = a, b = b, seed = d$seed),
    REF$gi.natural,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(
      Q = d$Q, mu = d$mu, a = a + d$mu, b = b + d$mu,
      seed = d$seed
    ),
    REF$gi.mu,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(Q = d$Q, a = a, b = b, seed = d$seed, max.size = 40),
    REF$gi.maxsize,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(Q = d$Q, a = a, b = b, seed = d$seed, lim = 0.5),
    REF$gi.lim,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(Q = d$Q, a = a, b = b, seed = d$seed, ind = lind),
    REF$gi.ind,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(
      Q = d$Q, a = a, b = b, seed = d$seed, ind = lind,
      use.reordering = "limits"
    ),
    REF$gi.limits,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(
      Q = d$Q, a = a, b = b, seed = d$seed,
      use.reordering = "sparsity"
    ),
    REF$gi.sparsity,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(Q.chol = chol(d$Q), a = a, b = b, seed = d$seed),
    REF$gi.chol,
    tolerance = tol
  )
  expect_equal(gaussint.ref(Q = d$Q, a = a, b = b, seed = 7),
    REF$gi.seed1,
    tolerance = tol
  )
  expect_equal(
    gaussint.ref(Q = d$Q, a = a, b = rep(Inf, d$n), seed = d$seed),
    REF$gi.onesided,
    tolerance = tol
  )
  ## Very narrow intervals, where samples can get negligible weights
  sd <- sqrt(excursions.variances(Q = d$Q))
  k <- seq(5, 95, by = 10)
  a.n <- a
  b.n <- b
  a.n[k] <- d$mu[k] + 0.3 * sd[k]
  b.n[k] <- a.n[k] + 1e-4 * sd[k]
  res <- gaussint.ref(Q = d$Q, a = a.n, b = b.n, seed = d$seed)
  expect_equal(res, REF$gi.narrow, tolerance = tol)
  ## The probability is tiny, and expect_equal() uses an absolute tolerance
  ## for values smaller than the tolerance, so compare on the log scale
  expect_equal(log(res[c("P", "E")]), log(REF$gi.narrow[c("P", "E")]),
    tolerance = tol
  )
  expect_equal(excursions:::excursions.rand(5, seed = 1:6, n.threads = 1),
    REF$rand,
    tolerance = tol
  )
})

test_that("Integration is exact for independent components", {
  ## With a diagonal precision matrix all samples get the same weight
  n <- 8
  sd <- seq(0.5, 2, length.out = n)
  Q <- Matrix::Diagonal(n, 1 / sd^2)
  a <- seq(-2, -1, length.out = n)
  b <- seq(0.5, 3, length.out = n)
  res <- gaussint(Q = Q, a = a, b = b, seed = 1:6, max.threads = 1, n.iter = 100)
  p <- pnorm(b / sd) - pnorm(a / sd)
  expect_equal(res$Pv, rev(cumprod(rev(p))), tolerance = 1e-12)
  expect_equal(res$P, prod(p), tolerance = 1e-12)
  ## Zero up to rounding errors
  expect_true(all(res$Ev < 1e-8))
})

test_that("Integration with integer and logical indices agree", {
  d <- testdata.spde(10)
  a <- d$mu - 2.5
  b <- d$mu + 2.5
  ind <- 30:70
  lind <- seq_len(d$n) %in% ind
  r.int <- gaussint(Q = d$Q, a = a, b = b, ind = ind, seed = d$seed, max.threads = 1, n.iter = 1000)
  r.log <- gaussint(Q = d$Q, a = a, b = b, ind = lind, seed = d$seed, max.threads = 1, n.iter = 1000)
  a.inf <- a
  b.inf <- b
  a.inf[!lind] <- -Inf
  b.inf[!lind] <- Inf
  r.inf <- gaussint(Q = d$Q, a = a.inf, b = b.inf, seed = d$seed, max.threads = 1, n.iter = 1000)
  expect_identical(r.int, r.log)
  expect_identical(r.int, r.inf)
  ## Leaving out nodes increases the probability
  r.all <- gaussint(Q = d$Q, a = a, b = b, seed = d$seed, max.threads = 1, n.iter = 1000)
  expect_gt(r.int$P, r.all$P)
})

test_that("Integration is reproducible", {
  d <- testdata.spde(10)
  a <- d$mu - 2.5
  b <- d$mu + 2.5
  run <- function(threads, seed = d$seed) {
    gaussint(Q = d$Q, a = a, b = b, seed = seed, max.threads = threads, n.iter = 2000)
  }
  r1 <- run(1)
  expect_identical(run(1), r1)
  ## With several threads, each thread uses its own random stream and a fixed
  ## block of samples, so the results only depend on the number of threads
  r2 <- run(2)
  expect_identical(run(2), r2)
  expect_equal(r2$P, r1$P, tolerance = 5 * r1$E / r1$P)
  ## A different seed gives different samples
  expect_false(identical(run(1, seed = 2:7)$Pv, r1$Pv))
})

test_that("Integration with a Cholesky factor matches the precision matrix", {
  d <- testdata.spde(10)
  a <- d$mu - 2.5
  b <- d$mu + 2.5
  r.Q <- gaussint(Q = d$Q, a = a, b = b, seed = d$seed, max.threads = 1, n.iter = 1000)
  r.L <- gaussint(Q.chol = chol(d$Q), a = a, b = b, seed = d$seed, max.threads = 1, n.iter = 1000)
  expect_equal(r.L, r.Q, tolerance = 1e-12)
})
