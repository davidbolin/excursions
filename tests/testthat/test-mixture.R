## Quantile of a Gaussian mixture at one location, computed with uniroot
mixture.quantile.ref <- function(p, mu, sd, w, br = c(-1000, 1000)) {
  stats::uniroot(
    function(x) sum(w * pnorm(x, mu, sd)) - p, br,
    tol = 1e-13
  )$root
}

test_that("Mixture quantiles are accurate", {
  set.seed(1)
  for (K in c(1, 2, 5, 20)) {
    n <- 200
    mu <- matrix(rnorm(K * n, sd = 3), K)
    sd <- matrix(exp(rnorm(K * n)), K)
    w <- runif(K)
    w <- w / sum(w)
    for (p in c(1e-6, 0.025, 0.5, 0.975, 1 - 1e-6)) {
      q <- excursions:::Fmix_inv_vec(p, mu, sd, w)
      F.q <- colSums(w * pnorm((rep(q, each = K) - mu) / sd))
      expect_equal(F.q, rep(p, n), tolerance = 1e-8)
      q.ref <- vapply(seq_len(20), function(i) {
        mixture.quantile.ref(p, mu[, i], sd[, i], w)
      }, 0.0)
      expect_equal(q[1:20], q.ref, tolerance = 1e-9)
    }
  }
})

test_that("Mixture quantiles are accurate for multimodal mixtures", {
  ## Newton's method oscillates for this mixture without a safeguard
  mu <- matrix(c(4.7230733, 2.2279823, 1.4794166, -4.4387700, -4.00706147))
  sd <- matrix(c(0.3021312, 0.1284623, 2.1579211, 0.9652709, 1.83982682))
  w <- c(0.2589725, 0.1430795, 0.4019409, 0.1387083, 0.05729875)
  q <- excursions:::Fmix_inv_vec(0.5, mu, sd, w)
  expect_equal(q, mixture.quantile.ref(0.5, mu, sd, w), tolerance = 1e-9)
  ## Converged Newton steps must not be replaced by bisection steps
  mu <- matrix(c(0.08693678, 0.28693678, -0.01306322))
  sd <- matrix(c(0.3791568, 0.3325425, 0.4531790))
  w <- c(0.5, 0.3, 0.2)
  for (p in c(0.05, 0.95)) {
    q <- excursions:::Fmix_inv_vec(p, mu, sd, w)
    expect_equal(q, mixture.quantile.ref(p, mu, sd, w), tolerance = 1e-9)
  }
})

test_that("Mixture quantiles in special cases", {
  ## Identical components
  q <- excursions:::Fmix_inv_vec(0.1, matrix(1, 3, 2), matrix(2, 3, 2), c(0.2, 0.3, 0.5))
  expect_equal(q, rep(1 + 2 * qnorm(0.1), 2))
  ## Quantiles outside the bounds are set to the bounds
  mu <- matrix(c(0, 5, -5), 1)
  sd <- matrix(c(1, 100, 100), 1)
  br <- c(-50, 50)
  q <- excursions:::Fmix_inv_vec(0.01, mu, sd, 1, br = br)
  expect_equal(q[1], qnorm(0.01))
  expect_identical(q[2:3], c(br[1], br[1]))
  q <- excursions:::Fmix_inv_vec(0.99, mu, sd, 1, br = br)
  expect_identical(q[2:3], c(br[2], br[2]))
})

test_that("Mixture confidence bands with sampling", {
  d <- testdata.mixture()
  run <- function(...) {
    simconf.mixture(
      alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w, n.iter = 2000,
      max.threads = 1, seed = 1, ...
    )
  }
  r <- run()
  expect_identical(run(), r)
  check.mixture.band(r, d)

  ## The marginal bounds are the quantiles of the mixture
  vars <- lapply(d$Q, function(Q) diag(solve(Q)))
  i <- 7
  expect_equal(r$a.marginal[i], mixture.quantile.ref(
    0.05, vapply(d$mu, `[`, 0, i), sqrt(vapply(vars, `[`, 0, i)), d$w
  ), tolerance = 1e-8)

  ## Subsets of the nodes
  ind <- 10:30
  r.ind <- run(ind = ind)
  expect_length(r.ind$a, length(ind))
  expect_equal(r.ind$a.marginal, r$a.marginal[ind])
  expect_true(all(r.ind$b - r.ind$a <= r$b[ind] - r$a[ind]))

  ## A single component
  r1 <- simconf.mixture(
    alpha = 0.1, mu = d$mu[1], Q = d$Q[1], w = 1, n.iter = 2000,
    max.threads = 1, seed = 1
  )
  expect_equal(r1$a.marginal, d$mu[[1]] + qnorm(0.05) * sqrt(vars[[1]]),
    tolerance = 1e-8
  )
  expect_true(all(r1$a < r1$a.marginal))
})

test_that("Mixture confidence bands with sequential integration", {
  d <- testdata.mixture()
  run <- function(...) {
    simconf.mixture(
      alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w, n.iter = 2000,
      max.threads = 1, seed = d$seed, mix.samp = FALSE, ...
    )
  }
  ## Used to fail when ind was not given
  r <- run()
  expect_identical(run(), r)
  check.mixture.band(r, d)

  ## Subsets of the nodes, where the nodes are reordered internally.
  ## The variances used to be permuted twice in this case.
  ind <- 10:30
  vars <- lapply(d$Q, function(Q) diag(solve(Q)))
  r.ind <- run(ind = ind)
  r.vars <- run(ind = ind, vars = vars)
  expect_equal(r.ind$a.marginal, r$a.marginal[ind], tolerance = 1e-10)
  expect_equal(r.ind$b.marginal, r$b.marginal[ind], tolerance = 1e-10)
  expect_equal(r.vars$a.marginal, r.ind$a.marginal, tolerance = 1e-10)
  expect_equal(r.vars$a, r.ind$a, tolerance = 1e-10)
  expect_equal(r.ind$vars, r$vars[ind], tolerance = 1e-10)
  expect_identical(run(ind = seq_len(d$n) %in% ind), r.ind)
})

test_that("Mixture confidence bands check their input", {
  d <- testdata.mixture()
  expect_error(
    simconf.mixture(alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w[1:2]),
    "different length"
  )
  expect_error(
    simconf.mixture(
      alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w,
      vars = list(1, 2)
    ),
    "different length"
  )
})

test_that("Mixture confidence bands are unchanged", {
  ## Computed after the quantiles were computed to full precision, which
  ## changed the results by about 3e-5 compared to earlier versions
  d <- testdata.mixture()
  r <- simconf.mixture(
    alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w, n.iter = 2000,
    max.threads = 1, seed = d$seed, mix.samp = FALSE
  )
  expect_equal(
    c(
      a = ref.summary(r$a), b = ref.summary(r$b),
      am = ref.summary(r$a.marginal), bm = ref.summary(r$b.marginal)
    ),
    REF.MIXTURE$int,
    tolerance = 1e-8
  )
  r <- simconf.mixture(
    alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w, n.iter = 2000,
    max.threads = 1, seed = 1
  )
  expect_equal(
    c(
      a = ref.summary(r$a), b = ref.summary(r$b),
      am = ref.summary(r$a.marginal), bm = ref.summary(r$b.marginal)
    ),
    REF.MIXTURE$samp,
    tolerance = 1e-8
  )
})

test_that("simconf.mixture with an adaptive number of iterations", {
  d <- testdata.mixture()
  args <- list(
    alpha = 0.1, mu = d$mu, Q = d$Q, w = d$w, max.threads = 1,
    seed = d$seed, mix.samp = FALSE
  )
  r0 <- do.call(simconf.mixture, c(args, list(n.iter = 1000)))
  r1 <- do.call(simconf.mixture, c(args, list(n.iter = 10000, tol = 1)))
  expect_identical(r1$a, r0$a)
  expect_identical(r1$b, r0$b)
  expect_equal(r1$meta$tol, 1)
})
