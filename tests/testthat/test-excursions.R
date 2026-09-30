test_that("Excursions, alpha = 1, type = >", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1
  )
  r <- c(
    2.463175e-15, 1.030394e-09, 7.734328e-06, 0.002534894, 0.07579555,
    0.4188056, 0.8192485, 0.9746894, 0.9984611, 0.9999625, 0.9999997
  )
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 1, type = <", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = "<",
    seed = data$seed, max.threads = 1
  )
  r <- c(
    0.9999997, 0.9999619, 0.9984783, 0.9746628, 0.819054,
    0.4196815, 0.07603522, 0.002548885, 7.762954e-06, 1.033062e-09, 2.47039e-15
  )
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 1, type = =", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "=",
    seed = data$seed, max.threads = 1
  )
  r <- c(
    7.381175e-07, 8.137438e-05, 0.003200957, 0.05128127, 0.3331441,
    0.6420603, 0.1815013, 0.02201733, 0.001153959, 2.520208e-05,
    1.945911e-07
  )

  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 1, type = !=", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "!=",
    seed = data$seed, max.threads = 1
  )
  r <- c(
    0.9999993, 0.9999186, 0.996799, 0.9487187, 0.6668559, 0.3579397,
    0.8184987, 0.9779827, 0.998846, 0.9999748, 0.9999998
  )

  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = >", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1
  )
  r <- c(0, 0, 0, 0, 0, 0, 0, 0.9801319, 0.9988921, 0.9999753, 0.9999998)
  res$F[is.na(res$F)] <- 0
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = <", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "<",
    seed = data$seed, max.threads = 1
  )
  ## Before the two nodes with the largest marginal probabilities were
  ## ordered by their probabilities, F[1] was smaller than F[2]
  r <- c(0.9999995, 0.9999425, 0.9979041, 0.9679957, 0, 0, 0, 0, 0, 0, 0)
  res$F[is.na(res$F)] <- 0
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = =", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "=",
    seed = data$seed, max.threads = 1
  )
  r <- c(
    7.381175e-07, 8.137438e-05, 0.003200957, 0.05128127, 1, 1,
    1, 0.02201733, 0.001153959, 2.520208e-05, 1.945911e-07
  )
  res$F[is.na(res$F)] <- 1
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = !=", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "!=",
    seed = data$seed, max.threads = 1
  )
  r <- c(
    0.9999993, 0.9999186, 0.996799, 0.9487187, 0, 0, 0,
    0.9779827, 0.998846, 0.9999748, 0.9999998
  )
  res$F[is.na(res$F)] <- 0
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, move u to mu", {
  data <- integration.testdata1()

  res <- excursions(
    alpha = 0.1, u = 1, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1
  )
  res2 <- excursions(
    alpha = 0.1, u = 0, mu = data$mu - 1, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1
  )
  res$F[is.na(res$F)] <- 0
  res2$F[is.na(res2$F)] <- 0
  expect_equal(res$F, res2$F, tolerance = 1e-7)
})

test_that("Excursions, input variances", {
  data <- integration.testdata1()

  vars <- diag(solve(data$Q))
  res1 <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, vars = vars, max.threads = 1
  )
  res2 <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1
  )
  expect_equal(res1$F, res2$F, tolerance = 1e-7)
})

test_that("Excursions, ind argument order", {
  data <- integration.testdata1()

  vars <- diag(solve(data$Q))

  ind1 <- c(1, 2, 3, 4)
  ind2 <- c(4, 3, 2, 1)
  ind3 <- rep(FALSE, length(data$mu))
  ind3[1:4] <- TRUE
  res1 <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, ind = ind1, max.threads = 1
  )
  res2 <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, ind = ind2, max.threads = 1
  )
  res3 <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, ind = ind3, max.threads = 1
  )

  expect_equal(res1$F, res2$F, tolerance = 1e-7)
  expect_equal(res2$F, res3$F, tolerance = 1e-7)
})

# Tests to add:

# Test that Q.chol and Q gives the same result

# Test the max.size argument

# Test reo

# Test rho

# Test max.threads

# Test QC method

# Test ind argument


# res1 = excursions2::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', seed = seed, max.threads = 1)
# res2 = excursions::excursions(alpha=0.1, u=1, mu=mu.x+0.1, Q=Q.x, type='!=', max.threads = 1)

# plot(res1$F)
# lines(res2$F, col=2)


test_that("Excursion sets are unchanged", {
  d <- testdata.spde(10)
  tol <- 1e-8
  for (type in c(">", "<", "=", "!=")) {
    r <- excursions(
      alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = type,
      seed = d$seed, max.threads = 1, n.iter = 2000
    )
    expect_equal(
      c(
        F = ref.summary(r$F), E = ref.summary(r$E), M = ref.summary(r$M),
        vars = ref.summary(r$vars)
      ),
      REF[[paste0("ex", type)]],
      tolerance = tol
    )
  }
  r <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = ">",
    seed = d$seed, max.threads = 1, n.iter = 2000, ind = 20:80
  )
  expect_equal(c(F = ref.summary(r$F), E = ref.summary(r$E)), REF$ex.ind,
    tolerance = tol
  )
  rho <- pnorm(d$mu - 0.5, sd = 1.2 * sqrt(excursions.variances(Q = d$Q)))
  r <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = ">",
    seed = d$seed, max.threads = 1, n.iter = 2000, method = "QC", rho = rho
  )
  expect_equal(c(F = ref.summary(r$F), E = ref.summary(r$E)), REF$ex.qc,
    tolerance = tol
  )
})

test_that("Excursions with a Cholesky factor matches the precision matrix", {
  ## The Cholesky factor used to be treated as a precision matrix when the
  ## nodes were reordered. The reordering is computed from the sparsity
  ## pattern of the matrix that is given, so use the same one for both.
  d <- testdata.spde(8)
  for (type in c(">", "=")) {
    r.Q <- excursions(
      alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = type,
      seed = d$seed, max.threads = 1, n.iter = 1000
    )
    expect_false(identical(r.Q$meta$reo, seq_len(d$n)))
    r.L <- excursions(
      alpha = 0.1, u = 0.5, mu = d$mu, Q.chol = chol(d$Q), type = type,
      seed = d$seed, max.threads = 1, n.iter = 1000, reo = r.Q$meta$reo
    )
    expect_equal(r.L$F, r.Q$F, tolerance = 1e-10)
    expect_equal(r.L$vars, r.Q$vars, tolerance = 1e-10)
  }
  ## With the natural ordering, the factor is used directly
  r.Q <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = ">",
    seed = d$seed, max.threads = 1, n.iter = 1000, reo = seq_len(d$n)
  )
  r.L <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q.chol = chol(d$Q), type = ">",
    seed = d$seed, max.threads = 1, n.iter = 1000, reo = seq_len(d$n)
  )
  expect_equal(r.L$F, r.Q$F, tolerance = 1e-10)
  ## With the computed ordering, the results agree up to Monte Carlo error
  r.Q <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = ">",
    seed = d$seed, max.threads = 1, n.iter = 1000
  )
  r.L <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q.chol = chol(d$Q), type = ">",
    seed = d$seed, max.threads = 1, n.iter = 1000
  )
  ok <- !is.na(r.Q$F) & !is.na(r.L$F)
  expect_true(all(abs(r.L$F[ok] - r.Q$F[ok]) < 0.01))
})

test_that("Nodes above the threshold are integrated in order of their probabilities", {
  d <- testdata.spde(10)
  set.seed(3)
  rho <- runif(d$n)
  for (alpha in c(0.05, 0.5, 1)) {
    reo <- excursions:::excursions.permutation(rho, NULL, TRUE, alpha, d$Q)
    above <- rho > 1 - alpha
    ## The nodes above the threshold come last, in increasing order of rho
    expect_setequal(tail(reo, sum(above)), which(above))
    expect_equal(tail(reo, sum(above)), which(above)[order(rho[above])])
    expect_setequal(reo, seq_len(d$n))
  }
  ## All nodes above the threshold
  expect_equal(
    excursions:::excursions.permutation(rho, NULL, TRUE, 1.5, d$Q),
    order(rho)
  )
})

test_that("Constrained reordering with merged constraint sets", {
  ## Merging runs of single node constraint sets must give the same ordering
  ## as calling CAMD with all the sets
  camd <- function(Q, cind) {
    Q <- excursions:::private.as.dgCMatrix(Q)
    out <- .C("reordering",
      nin = as.integer(nrow(Q)), Mp = as.integer(Q@p),
      Mi = as.integer(Q@i), reo = integer(nrow(Q)),
      cind = as.integer(cind), PACKAGE = "excursions"
    )
    out$reo + 1L
  }
  d <- testdata.spde(10)
  n <- d$n
  set.seed(4)
  for (rep in 1:5) {
    ## A large set, single node sets, a set with a few nodes, and more
    ## single node sets
    p <- sample(n)
    cind <- integer(n)
    cind[p[61:80]] <- 1:20
    cind[p[81:85]] <- 21
    cind[p[86:n]] <- 22:(21 + n - 85)
    expect_identical(excursions:::private.camd(d$Q, cind), camd(d$Q, cind))
  }
})
