test_that("Monte Carlo results are unchanged", {
  d <- testdata.spde(10)
  X <- testdata.spde.samples(d)
  tol <- 1e-10

  r <- excursions.mc(X, u = 0.5, type = ">", alpha = 0.1)
  expect_equal(
    c(F = ref.summary(r$F), E = ref.summary(r$E), M = ref.summary(r$M)),
    REF$exmc,
    tolerance = tol
  )
  r <- excursions.mc(X, u = 0.5, type = "=", alpha = 0.1, ind = 20:80)
  expect_equal(
    c(F = ref.summary(r$F), E = ref.summary(r$E), M = ref.summary(r$M)),
    REF$exmc.ind,
    tolerance = tol
  )
  r <- contourmap.mc(X,
    n.levels = 3, alpha = 0.1,
    compute = list(F = TRUE, measures = c("P1", "P2"))
  )
  expect_equal(
    c(
      F = ref.summary(r$F), E = ref.summary(r$E), M = ref.summary(r$M),
      P1 = r$P1, P2 = r$P2
    ),
    REF$cmmc,
    tolerance = tol
  )
  r <- simconf.mc(X, alpha = 0.1)
  expect_equal(
    c(
      a = ref.summary(r$a), b = ref.summary(r$b),
      am = ref.summary(r$a.marginal), bm = ref.summary(r$b.marginal)
    ),
    REF$scmc,
    tolerance = tol
  )
  r <- simconf.mc(X, alpha = 0.1, ind = 20:80)
  expect_equal(c(a = ref.summary(r$a), b = ref.summary(r$b)), REF$scmc.ind,
    tolerance = tol
  )
})

test_that("simconf.mc prints nothing", {
  d <- testdata.spde(6)
  X <- testdata.spde.samples(d, 100)
  expect_silent(simconf.mc(X, alpha = 0.1))
})

test_that("mcint agrees with the direct computation", {
  d <- testdata.spde(6)
  X <- testdata.spde.samples(d, 300)
  a <- d$mu - 1.5
  b <- d$mu + 1.5
  inside <- a < X & X < b
  ## P(a_j < X_j < b_j for all j >= i)
  P.ref <- vapply(seq_len(d$n), function(i) {
    mean(colSums(!inside[i:d$n, , drop = FALSE]) == 0)
  }, 0.0)
  r <- excursions:::mcint(X = X, a = a, b = b)
  expect_equal(r$Pv, P.ref)
  expect_equal(r$P, P.ref[1])

  ## Integer and logical indices
  ind <- 5:20
  lind <- seq_len(d$n) %in% ind
  r.int <- excursions:::mcint(X = X, a = a, b = b, ind = ind)
  r.log <- excursions:::mcint(X = X, a = a, b = b, ind = lind)
  expect_identical(r.int, r.log)
  P.ind <- vapply(seq_len(d$n), function(i) {
    mean(colSums(!inside[intersect(i:d$n, ind), , drop = FALSE]) == 0)
  }, 0.0)
  expect_equal(r.int$Pv, P.ind)

  ## Missing values give missing probabilities from that node down
  X.na <- X
  X.na[10, 1] <- NA
  r.na <- excursions:::mcint(X = X.na, a = a, b = b)
  expect_true(all(is.na(r.na$Pv[1:10])))
  expect_equal(r.na$Pv[11:d$n], P.ref[11:d$n])
})

test_that("Row quantiles agree with quantile()", {
  set.seed(1)
  X <- matrix(rnorm(60), 6, 10)
  X[2, ] <- round(X[2, ]) ## ties
  X[3, ] <- 1 ## constant
  sorted <- excursions:::private.sort.rows(X)
  for (p in c(0, 0.01, 0.05, 1 / 9, 0.5, 0.95, 1)) {
    expect_identical(
      excursions:::private.row.quantile(sorted, p),
      unname(apply(X, 1, quantile, probs = p))
    )
  }
  ## One and two samples
  for (m in 1:2) {
    Xm <- X[, seq_len(m), drop = FALSE]
    sm <- excursions:::private.sort.rows(Xm)
    expect_equal(dim(sm), dim(Xm))
    for (p in c(0, 0.3, 1)) {
      expect_identical(
        excursions:::private.row.quantile(sm, p),
        unname(apply(Xm, 1, quantile, probs = p))
      )
    }
  }
  ## Missing values are not sorted
  X[1, 1] <- NA
  expect_null(excursions:::private.sort.rows(X))
})

test_that("The coverage for simconf.mc agrees with and without sorting", {
  d <- testdata.spde(6)
  X <- testdata.spde.samples(d, 200)
  sorted <- excursions:::private.sort.rows(X)
  for (x in c(0.001, 0.01, 0.05)) {
    expect_identical(
      excursions:::fsamp.opt(x, X, sorted = sorted),
      excursions:::fsamp.opt(x, X)
    )
  }
})
