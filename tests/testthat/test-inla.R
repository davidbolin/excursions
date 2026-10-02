test_that("The mode configuration is the first one with the largest posterior", {
  result <- list(misc = list(configs = list(config = lapply(
    c(-3, -1, -2, -1),
    function(lp) list(log.posterior = lp)
  ))))
  expect_equal(excursions:::private.mode.config(result), 2)
})

## Keep only the k configurations with the largest posterior, to reduce the
## number of refits for iNIQC
testdata.inla.trim <- function(result, k = 3) {
  configs <- result$misc$configs$config
  lp <- vapply(configs, function(x) x$log.posterior, 0.0)
  keep <- sort(order(lp, decreasing = TRUE)[seq_len(k)])
  result$misc$configs$config <- configs[keep]
  result$misc$configs$nconfig <- k
  result
}

test_that("Marginal probabilities of intervals", {
  skip_on_cran()
  local_exc_safe_inla()
  data <- testdata.inla()
  ## The same as inla.pmarginal, up to its error from the numerical
  ## integration with 2048 points
  i <- c(1, 5)
  ref <- function(margs) {
    t(vapply(margs, function(m) {
      c(INLA::inla.pmarginal(-1, m), INLA::inla.pmarginal(2, m))
    }, numeric(2)))
  }
  p <- excursions:::inla.get.marginal.int(i, a = -1, b = 2, result = data$result)
  expect_equal(dim(p), c(2, 2))
  expect_equal(p, ref(data$result$marginals.linear.predictor[i]),
    tolerance = 2e-3, ignore_attr = TRUE
  )
  p <- excursions:::inla.get.marginal.int(i,
    a = -1, b = 2, result = data$result,
    effect.name = "ar"
  )
  expect_equal(p, ref(data$result$marginals.random$ar[i]),
    tolerance = 2e-3, ignore_attr = TRUE
  )
})

test_that("contourmap.inla with the QC method", {
  skip_on_cran()
  local_exc_safe_inla()
  data <- testdata.inla()
  run <- function(method) {
    contourmap.inla(data$result, data$stack,
      tag = "pred", method = method, n.levels = 2, alpha = 0.1,
      seed = data$seed, max.threads = 1, compute = list(F = TRUE)
    )
  }
  r.qc <- run("QC")
  r.eb <- run("EB")
  expect_true(all(r.qc$F >= 0 & r.qc$F <= 1, na.rm = TRUE))
  ## The likelihood is Gaussian, so the methods give similar results
  expect_equal(r.qc$F, r.eb$F, tolerance = 0.05)
})

test_that("excursions.inla with the iNIQC method refits the model", {
  skip_on_cran()
  local_exc_safe_inla()
  data <- testdata.inla()
  result <- testdata.inla.trim(data$result)

  ## Record whether the marginals come from the original fit or a refit
  from.refit <- logical(0)
  get.marginal <- excursions:::inla.get.marginal
  local_mocked_bindings(inla.get.marginal = function(i, u, result, ...) {
    from.refit <<- c(
      from.refit,
      rep(isTRUE(result$.args$control.mode$fixed), length(i))
    )
    get.marginal(i, u = u, result = result, ...)
  })
  r.in <- excursions.inla(result, data$stack,
    tag = "pred", method = "iNIQC",
    u = 0, type = ">", seed = data$seed, max.threads = 1
  )
  n.ind <- length(INLA::inla.stack.index(data$stack, "pred")$data)
  ## The marginals of each node from the original fit, and from the refit
  ## of each configuration
  expect_equal(sum(!from.refit), n.ind)
  expect_equal(sum(from.refit), 3 * n.ind)

  expect_true(all(r.in$F >= 0 & r.in$F <= 1))
  r.ni <- excursions.inla(result, data$stack,
    tag = "pred", method = "NIQC",
    u = 0, type = ">", seed = data$seed, max.threads = 1
  )
  expect_equal(r.in$F, r.ni$F, tolerance = 0.05)
})

test_that("excursions.inla with the iNIQC method for a random effect", {
  skip_on_cran()
  local_exc_safe_inla()
  data <- testdata.inla()
  result <- testdata.inla.trim(data$result)
  r.in <- excursions.inla(result,
    name = "ar", method = "iNIQC",
    u = 0, type = ">", seed = data$seed, max.threads = 1
  )
  r.ni <- excursions.inla(result,
    name = "ar", method = "NIQC",
    u = 0, type = ">", seed = data$seed, max.threads = 1
  )
  expect_true(all(r.in$F >= 0 & r.in$F <= 1))
  expect_equal(r.in$F, r.ni$F, tolerance = 0.05)
})

test_that("simconf.inla", {
  skip_on_cran()
  local_exc_safe_inla()
  data <- testdata.inla()
  run <- function(method, ...) {
    simconf.inla(data$result, data$stack,
      tag = "pred", method = method, alpha = 0.1,
      seed = data$seed, max.threads = 1, ...
    )
  }
  check <- function(r) {
    expect_true(all(r$a < r$a.marginal))
    expect_true(all(r$a.marginal < r$b.marginal))
    expect_true(all(r$b.marginal < r$b))
  }
  r.eb <- run("EB")
  check(r.eb)
  r.ni <- run("NI")
  check(r.ni)
  r.ni.int <- run("NI", inla.sample = FALSE)
  check(r.ni.int)
  expect_equal(r.ni.int$a, r.ni$a, tolerance = 0.01)
  expect_equal(r.ni.int$b, r.ni$b, tolerance = 0.01)
  expect_equal(r.ni.int$a.marginal, r.ni$a.marginal, tolerance = 1e-10)

  ## Integer indices are used for the band, which is narrower for fewer nodes
  ind <- 3:6
  r.ind <- run("EB", ind = ind)
  expect_length(r.ind$a, length(ind))
  expect_true(all(r.ind$b - r.ind$a < r.eb$b[ind] - r.eb$a[ind]))
  expect_equal(r.ind$a.marginal, r.eb$a.marginal[ind], tolerance = 1e-10)
})

test_that("Distribution functions of marginals", {
  ## Gaussian densities on grids like those of INLA, where the distribution
  ## function, normalised on the range of the grid, is known
  set.seed(1)
  mu <- rnorm(20)
  sd <- exp(rnorm(20, -1, 0.5))
  margs <- lapply(seq_along(mu), function(k) {
    x <- mu[k] + sd[k] * sort(c(qnorm(seq(0.001, 0.999, length.out = 40)), rnorm(3)))
    cbind(x = x, y = dnorm(x, mu[k], sd[k]))
  })
  q <- mu + sd * rnorm(20)
  lo <- vapply(margs, function(m) min(m[, 1]), 0)
  hi <- vapply(margs, function(m) max(m[, 1]), 0)
  exact <- (pnorm(pmin(pmax(q, lo), hi), mu, sd) - pnorm(lo, mu, sd)) /
    (pnorm(hi, mu, sd) - pnorm(lo, mu, sd))
  p <- excursions:::private.pmarginals(margs, q)
  expect_equal(p, exact, tolerance = 1e-8)
  ## Marginals as lists, and of different lengths
  margs.l <- lapply(margs, function(m) list(x = m[, 1], y = m[, 2]))
  margs.l[[3]] <- lapply(margs.l[[3]], function(v) v[-(1:5)])
  p.l <- excursions:::private.pmarginals(margs.l, q)
  expect_equal(p.l[-3], p[-3])
  expect_equal(
    p.l[3],
    excursions:::private.pmarginals(margs.l[3], q[3])
  )
  ## q outside the range, and recycled
  expect_equal(excursions:::private.pmarginals(margs, -100), rep(0, 20))
  expect_equal(excursions:::private.pmarginals(margs, 100), rep(1, 20))
  ## Points with negligible or missing densities are removed
  m <- margs[[1]]
  m2 <- rbind(m, c(max(m[, 1]) + 1, 0), c(max(m[, 1]) + 2, NA))
  expect_equal(
    excursions:::private.pmarginals(list(m2), q[1]),
    excursions:::private.pmarginals(list(m), q[1])
  )
})
