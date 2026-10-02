test_that("Excursions, alpha = 1, type = >", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1, n.iter = 10000
  )
  r <- c(
    2.453585944e-15, 1.028039605e-09, 7.741444498e-06, 0.002538569057,
    0.07582882438, 0.4188539218, 0.819121314, 0.9746501927, 0.9984674663,
    0.9999620956, 0.9999996732
  )
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 1, type = <", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = "<",
    seed = data$seed, max.threads = 1, n.iter = 10000
  )
  r <- c(
    0.9999996732, 0.9999622207, 0.9984746347, 0.9746765193, 0.8191915657,
    0.4197256656, 0.07604785403, 0.002541501738, 7.755477415e-06,
    1.029995247e-09, 2.461995094e-15
  )
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 1, type = =", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "=",
    seed = data$seed, max.threads = 1, n.iter = 10000
  )
  r <- c(
    7.38117649e-07, 8.185870815e-05, 0.003203209453, 0.05128514472,
    0.3325649633, 0.6421184302, 0.1810603749, 0.02201765987,
    0.001159227166, 2.54498799e-05, 1.945910227e-07
  )

  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 1, type = !=", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "!=",
    seed = data$seed, max.threads = 1, n.iter = 10000
  )
  r <- c(
    0.9999992619, 0.9999181413, 0.9967967905, 0.9487148553, 0.6674350367,
    0.3578815698, 0.8189396251, 0.9779823401, 0.9988407728, 0.9999745501,
    0.9999998054
  )

  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = >", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1, size.tol = NULL,
    n.iter = 10000
  )
  r <- c(
    0, 0, 0, 0, 0, 0, 0, 0.9800983735, 0.9988969148, 0.9999750936,
    0.9999998054
  )
  res$F[is.na(res$F)] <- 0
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = <", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "<",
    seed = data$seed, max.threads = 1, size.tol = NULL,
    n.iter = 10000
  )
  ## Before the two nodes with the largest marginal probabilities were
  ## ordered by their probabilities, F[1] was smaller than F[2]
  r <- c(
    0.9999994565, 0.9999430285, 0.9978992647, 0.9680120958, 0, 0, 0, 0,
    0, 0, 0
  )
  res$F[is.na(res$F)] <- 0
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = =", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "=",
    seed = data$seed, max.threads = 1, size.tol = NULL,
    n.iter = 10000
  )
  r <- c(
    7.38117649e-07, 8.185870815e-05, 0.003203209453, 0.05128514472, 1, 1,
    1, 0.02201765987, 0.001159227166, 2.54498799e-05, 1.945910227e-07
  )
  res$F[is.na(res$F)] <- 1
  expect_equal(res$F, r, tolerance = 1e-7)
})

test_that("Excursions, alpha = 0.1, type = !=", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 0.1, u = 0, mu = data$mu + 0.1, Q = data$Q, type = "!=",
    seed = data$seed, max.threads = 1, size.tol = NULL,
    n.iter = 10000
  )
  r <- c(
    0.9999992619, 0.9999181413, 0.9967967905, 0.9487148553, 0, 0, 0,
    0.9779823401, 0.9988407728, 0.9999745501, 0.9999998054
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
  ## With a fixed number of iterations
  d <- testdata.spde(10)
  tol <- 1e-8
  for (type in c(">", "<", "=", "!=")) {
    r <- excursions(
      alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = type,
      seed = d$seed, max.threads = 1, n.iter = 2000,
      size.tol = NULL
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
    seed = d$seed, max.threads = 1, n.iter = 2000,
      size.tol = NULL, ind = 20:80
  )
  expect_equal(c(F = ref.summary(r$F), E = ref.summary(r$E)), REF$ex.ind,
    tolerance = tol
  )
  rho <- pnorm(d$mu - 0.5, sd = 1.2 * sqrt(excursions.variances(Q = d$Q)))
  r <- excursions(
    alpha = 0.1, u = 0.5, mu = d$mu, Q = d$Q, type = ">",
    seed = d$seed, max.threads = 1, n.iter = 2000,
      size.tol = NULL, method = "QC", rho = rho
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

test_that("Excursions with adaptive number of iterations", {
  data <- integration.testdata1()
  res0 <- excursions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1, size.tol = NULL,
    n.iter = 10000
  )
  res1 <- excursions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1, tol = 1e-3
  )
  res2 <- excursions(
    alpha = 0.1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1, tol = 1e-4, n.iter = 1e5
  )
  expect_equal(res0$meta$n.iter.used, 10000)
  expect_equal(res1$meta$n.iter.used, 1000)
  expect_gt(res2$meta$n.iter.used, 1000)
  expect_lt(res2$meta$n.iter.used, 1e5)
  expect_equal(res1$E, res0$E)
  expect_equal(res2$E, res0$E)
  expect_equal(res1$F, res0$F, tolerance = 1e-3)
  ## The error is controlled at the boundary of the excursion set
  boundary <- which(res2$E == 1)[which.min(res2$F[res2$E == 1])]
  expect_lte(res2$meta$Fe[boundary], 1e-4)
})

test_that("Excursions with adaptive number of iterations and alpha = 1", {
  data <- integration.testdata1()
  res <- excursions(
    alpha = 1, u = 0, mu = data$mu, Q = data$Q, type = ">",
    seed = data$seed, max.threads = 1, tol = 2e-4, n.iter = 1e5
  )
  r <- c(
    2.463175e-15, 1.030394e-09, 7.734328e-06, 0.002534894, 0.07579555,
    0.4188056, 0.8192485, 0.9746894, 0.9984611, 0.9999625, 0.9999997
  )
  expect_gt(res$meta$n.iter.used, 1000)
  expect_lt(res$meta$n.iter.used, 1e5)
  expect_equal(res$F, r, tolerance = 2e-3)
  ## The error is controlled where F passes 0.5
  expect_lte(max(res$meta$Fe[6:7]), 2e-4)
})

test_that("Excursions with covariances fixes fewer nodes and gives the same results", {
  d <- testdata.spde(30)
  Q <- as(d$Q, "CsparseMatrix")
  sel <- excursions:::private.selected.inverse(Q)
  Qt <- as(Matrix::triu(Q), "TsparseMatrix")
  ## Upper triangle on the pattern of Q, as the Qinv of INLA
  Qinv <- Matrix::sparseMatrix(
    i = Qt@i + 1, j = Qt@j + 1,
    x = sel$cov(Qt@i + 1, Qt@j + 1), dims = dim(Q)
  )
  for (F.limit in c(0.1, 0.99)) {
    args <- list(
      alpha = 0.1, u = 0, mu = d$mu, Q = Q, type = ">", F.limit = F.limit,
      seed = d$seed, max.threads = 1, n.iter = 2000
    )
    ## All candidates in the order of rho
    r0 <- do.call(excursions, c(args, list(vars = sel$vars)))
    ## Covariances given, computed, or given as a symmetric matrix
    r1 <- do.call(excursions, c(args, list(vars = sel$vars, Qinv = Qinv)))
    r2 <- do.call(excursions, args)
    r3 <- do.call(excursions, c(args, list(Qinv = Matrix::forceSymmetric(Qinv, "U"))))
    for (r in list(r1, r2, r3)) {
      expect_equal(r$F, r0$F, tolerance = 1e-10)
      expect_identical(r$E, r0$E)
    }
    if (F.limit == 0.99) {
      ## Fewer nodes are in the order of rho
      expect_false(identical(r1$meta$reo, r0$meta$reo))
    }
  }
  ## The approximation of the reach is below the number of candidates
  marg <- excursions.marginals(type = ">", vars = sel$vars, mu = d$mu, u = 0)
  lims <- excursions.setlimits(marg, sel$vars, ">", FALSE, 0, d$mu)
  k <- excursions:::private.chain.reach(
    marg$rho, lims$a, lims$b, sel$vars, sel$cov, Q, 0.01
  )
  expect_lt(k, sum(marg$rho > 0.01) / 2)
  expect_gt(k, sum(!is.na(r0$F)) / 2)
})

test_that("Excursions adds nodes if the integration does not stop among the fixed nodes", {
  d <- testdata.spde(20)
  Q <- as(d$Q, "CsparseMatrix")
  vars <- excursions.variances(Q = Q)
  marg <- excursions.marginals(type = ">", vars = vars, mu = d$mu, u = 0)
  lims <- excursions.setlimits(marg, vars, ">", FALSE, 0, d$mu)
  run <- function(n.fixed) {
    excursions:::private.excursions.integrate(lims$a, lims$b, marg$rho, Q,
      is.chol = FALSE, F.limit = 0.99, m.size = d$n, n.fixed = n.fixed,
      n.iter = 2000, max.threads = 1, seed = d$seed, tol = NULL,
      tol.level = 0.9
    )
  }
  r0 <- run(Inf)
  r1 <- run(5)
  expect_gt(r1$n.fixed, 5)
  F0 <- F1 <- numeric(d$n)
  F0[r0$reo] <- r0$res$Pv
  F1[r1$reo] <- r1$res$Pv
  expect_equal(F1, F0, tolerance = 1e-10)
})

test_that("Excursions with Q.chol computes the covariances from the factor", {
  d <- testdata.spde(20)
  args <- list(
    alpha = 0.1, u = 0, mu = d$mu, type = ">", F.limit = 0.99,
    seed = d$seed, max.threads = 1, n.iter = 2000
  )
  r0 <- do.call(excursions, c(args, list(Q = d$Q)))
  r1 <- do.call(excursions, c(args, list(Q.chol = chol(d$Q))))
  expect_equal(r1$F, r0$F, tolerance = 1e-8)
  expect_identical(r1$E, r0$E)
})

test_that("Excursions checks Qinv", {
  d <- testdata.spde(10)
  expect_error(
    excursions(
      alpha = 0.1, u = 0, mu = d$mu, Q = d$Q, type = ">",
      Qinv = Matrix::Diagonal(5)
    ),
    "Qinv must be a sparse matrix"
  )
})

test_that("The reach estimate allows nodes without limits", {
  d <- testdata.spde(20)
  Q <- as(d$Q, "CsparseMatrix")
  sel <- excursions:::private.selected.inverse(Q)
  marg <- excursions.marginals(type = ">", vars = sel$vars, mu = d$mu, u = 0)
  lims <- excursions.setlimits(marg, sel$vars, ">", FALSE, 0, d$mu)
  k0 <- excursions:::private.chain.reach(
    marg$rho, lims$a, lims$b, sel$vars, sel$cov, Q, 0.01
  )
  ## Nodes with probability one have no limits, as in the QC method
  top <- order(marg$rho, decreasing = TRUE)[1:10]
  a <- lims$a
  a[top] <- -Inf
  k1 <- excursions:::private.chain.reach(
    marg$rho, a, lims$b, sel$vars, sel$cov, Q, 0.01
  )
  expect_false(is.na(k1))
  expect_gte(k1, k0)
  ## Two sided limits are not supported
  b <- lims$b
  b[top[1]] <- 10
  expect_true(is.na(excursions:::private.chain.reach(
    marg$rho, lims$a, b, sel$vars, sel$cov, Q, 0.01
  )))
})

test_that("Excursions with size.tol", {
  d <- testdata.spde(30)
  args <- list(alpha = 0.1, u = -1, mu = d$mu, Q = d$Q, type = ">", seed = d$seed)
  r.fixed <- do.call(excursions, c(args, list(size.tol = NULL, n.iter = 10000)))
  expect_equal(r.fixed$meta$n.iter.used, 10000)
  ## A large set needs fewer iterations for a target of 0.5%, and has almost
  ## the same size
  r <- do.call(excursions, c(args, list(size.tol = 0.005, n.iter = 10000)))
  expect_gt(sum(r$E), 200)
  expect_lt(r$meta$n.iter.used, 10000)
  expect_lte(abs(sum(r$E) - sum(r.fixed$E)), 0.02 * sum(r.fixed$E))
  ## The default is 0.1%, but at least half a node, with at most 20000
  ## iterations, so this set uses more iterations
  r <- do.call(excursions, args)
  expect_equal(r$meta$size.tol, 0.001)
  expect_gt(r$meta$n.iter.used, 10000)
  expect_lte(r$meta$n.iter.used, 20000)
  ## For a small set the target is half a node, and the size is close to
  ## the one with many iterations
  small <- list(alpha = 0.1, u = 1, mu = d$mu, Q = d$Q, type = ">", seed = d$seed)
  r <- do.call(excursions, c(small, list(n.iter = 3000, size.tol = 0.005)))
  r.many <- do.call(excursions, c(small, list(n.iter = 1e5, size.tol = NULL)))
  expect_gt(sum(r$E), 0)
  expect_lt(sum(r$E), 100)
  expect_lte(r$meta$n.iter.used, 3000)
  expect_lte(abs(sum(r$E) - sum(r.many$E)), 2)
  ## An empty set is exact after the first batch
  r <- excursions(
    alpha = 0.1, u = 10, mu = d$mu, Q = d$Q, type = ">", seed = d$seed
  )
  expect_equal(sum(r$E), 0)
  expect_equal(r$meta$n.iter.used, 1000)
  ## Without a set, n.iter iterations are used
  r <- excursions(
    alpha = 1, u = 0, mu = d$mu, Q = d$Q, type = ">", seed = d$seed,
    n.iter = 3000
  )
  expect_equal(r$meta$n.iter.used, 3000)
  expect_null(r$meta$size.tol)
})
