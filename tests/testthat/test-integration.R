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

## Run gaussint in a new R process with the given environment variables,
## such as "OMP_THREAD_LIMIT=2". The OpenMP runtime reads them when it starts,
## which on Linux is when R starts, so they are set in this process around the
## call and inherited by the new process. (The env argument of system2() is
## not supported for Rscript on Windows.)
gaussint.subprocess <- function(env, max.threads) {
  out <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(c(out, script)))
  writeLines(c(
    sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
    "suppressMessages(library(excursions))",
    "Q <- Matrix::forceSymmetric(Matrix::crossprod(Matrix::bandSparse(200,",
    "  k = 0:1, diagonals = list(rep(2, 200), rep(-1, 199)))))",
    "r <- gaussint(Q = Q, a = rep(-3, 200), b = rep(3, 200), seed = 1:6,",
    sprintf("  max.threads = %d, n.iter = 1001)", max.threads),
    "saveRDS(list(r = r, info = excursions:::private.openmp.info()),",
    sprintf("  %s)", deparse(out))
  ), script)
  vars <- character(0)
  if (length(env) > 0) {
    kv <- strsplit(env, "=", fixed = TRUE)
    vars <- stats::setNames(vapply(kv, `[`, "", 2), vapply(kv, `[`, "", 1))
  }
  status <- withr::with_envvar(vars, system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)),
    stdout = FALSE, stderr = FALSE
  ))
  if (status != 0 || !file.exists(out)) {
    return(NULL)
  }
  readRDS(out)
}

test_that("Integration respects the OpenMP thread limits", {
  skip_on_cran()
  skip_if_not(excursions:::private.openmp.info()[["openmp"]] == 1, "no OpenMP")
  ## More processors than the thread limit of 2 used below, so that asking for
  ## 4 threads asks for more threads than the limit allows
  skip_if(excursions:::private.openmp.info()[["num.procs"]] < 3, "fewer than 3 processors")

  ref2 <- gaussint.subprocess(character(0), max.threads = 2)
  skip_if(is.null(ref2), "could not run R in a subprocess")
  ref1 <- gaussint.subprocess(character(0), max.threads = 1)

  ## With OMP_THREAD_LIMIT=2, the runtime gives at most 2 threads. Asking for
  ## more used to leave some of the samples out of the estimate.
  lim <- gaussint.subprocess("OMP_THREAD_LIMIT=2", max.threads = 4)
  expect_equal(lim$info[["thread.limit"]], 2)
  expect_identical(lim$r, ref2$r)

  ## The default number of threads follows OMP_NUM_THREADS
  def <- gaussint.subprocess("OMP_NUM_THREADS=2", max.threads = 0)
  expect_identical(def$r, ref2$r)
  def1 <- gaussint.subprocess("OMP_NUM_THREADS=1", max.threads = 0)
  expect_identical(def1$r, ref1$r)
})

test_that("Integration with several threads agrees with one thread", {
  skip_if_not(excursions:::private.openmp.info()[["openmp"]] == 1, "no OpenMP")
  d <- testdata.spde(10)
  a <- d$mu - 2.5
  b <- d$mu + 2.5
  r1 <- gaussint(Q = d$Q, a = a, b = b, seed = d$seed, max.threads = 1, n.iter = 4000)
  r2 <- gaussint(Q = d$Q, a = a, b = b, seed = d$seed, max.threads = 2, n.iter = 4000)
  ## Different random streams, so equal up to the Monte Carlo error
  expect_false(identical(r1$Pv, r2$Pv))
  expect_lt(abs(r2$P - r1$P), 4 * sqrt(r1$E^2 + r2$E^2))
  ok <- r1$Pv > 0 & r2$Pv > 0
  expect_true(all(abs(r2$Pv[ok] - r1$Pv[ok]) < 4 * sqrt(r1$Ev[ok]^2 + r2$Ev[ok]^2) + 1e-12))
  expect_identical(
    gaussint(Q = d$Q, a = a, b = b, seed = d$seed, max.threads = 2, n.iter = 4000),
    r2
  )
})
