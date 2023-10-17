integration.testdata1 <- function() {
  n <- 11
  Q <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
  mu <- seq(-5, 5, length = n)
  L <- chol(Q)
  seed <- 1:6
  a <- rep(-3, n)
  b <- rep(3, n)
  list(n = n, Q = Q, mu = mu, seed = seed, L = L, a = a, b = b)
}


testdata.inla <- function(inla.mode = "classic") {
  if (require("INLA", quietly = TRUE)) {
    local_exc_safe_inla()
    set.seed(1)
    rho <- 0.99
    tau <- 15
    n <- 11
    n.obs <- 10
    x <- 1:n
    mu <- 10 * ((x < n / 2) * (x - n / 2) + (x >= n / 2) * (n / 2 - x) + n / 4) / n
    Q <- tau * sparseMatrix(
      i = c(1:n, 2:n), j = c(1:n, 1:(n - 1)),
      x = c(1, rep(1 + rho^2, n - 2), 1, rep(-rho, n - 1)),
      dims = c(n, n), symmetric = TRUE
    )

    X <- mu + INLA::inla.qsample(1, Q, seed = 12345L)[, 1]

    obs.loc <- c(3, 4, 6, 8, 2, 9, 5, 11, 7, 1)
    A <- sparseMatrix(i = 1:n.obs, j = obs.loc, x = rep(1, n.obs), dims = c(n.obs, n))
    Y <- as.vector(A %*% X +
      INLA::inla.qsample(1, Diagonal(n.obs, 1.0),
        seed = 98765L
      )[, 1])

    ef <- list(c(list(ar = x), list(cov = mu)))
    s.obs <- INLA::inla.stack(data = list(y = Y), A = list(A), effects = ef, tag = "obs")
    s.pre <- INLA::inla.stack(data = list(y = NA), A = list(1), effects = ef, tag = "pred")
    stack <- INLA::inla.stack(s.obs, s.pre)
    formula <- y ~ -1 + cov + f(ar, model = "ar1")
    result <- INLA::inla(
      formula = formula,
      data = INLA::inla.stack.data(stack),
      control.predictor = list(
        A = INLA::inla.stack.A(stack),
        compute = TRUE
      ),
      control.compute = list(
        config = TRUE,
        return.marginals.predictor = TRUE
      ),
      control.family = list(hyper = list(
        prec = list(initial = log(10), fixed = TRUE)
      )),
      num.threads = "1:1",
      inla.mode = inla.mode
    )

    seed <- 1:6
    return(list(result = result, stack = stack, seed = seed, n = n))
  }
}

testdata.inla.small <- function(inla.mode = "classic") {
  if (require("INLA", quietly = TRUE)) {
    local_exc_safe_inla()
    set.seed(1)
    rho <- 0.9
    sigma <- 3
    n <- 7
    n.obs <- 5
    x <- 1:n
    mu <- 10 * ((x < n / 2) * (x - n / 2) + (x >= n / 2) * (n / 2 - x) + n / 4) / n
    # Q <- sparseMatrix(i=c(1:n, 2:n), j=c(1:n, 1:(n-1)),
    #                      x=c(1,rep(1+rho^2, n-2),1, rep(-rho, n-1)),
    #                      dims=c(n, n), symmetric=TRUE) / (1 - rho^2) / sigma^2

    # X <- mu + INLA::inla.qsample(1, Q, seed=12345L)[,1]

    obs.loc <- c(3, 4, 6, 2, 5)
    A <- sparseMatrix(i = 1:n.obs, j = obs.loc, x = rep(1, n.obs), dims = c(n.obs, n))
    # Y <- as.vector(A %*% X +
    #                 INLA::inla.qsample(1, Diagonal(n.obs, 10.0),
    #                                    seed=98765L)[,1])
    Y <- c(4.812842, 5.165710, 2.786433, 1.391536, 3.094942)
    ef <- list(c(list(ar = x), list(cov = mu)))
    s.obs <- INLA::inla.stack(data = list(y = Y), A = list(A), effects = ef, tag = "obs")
    s.pre <- INLA::inla.stack(data = list(y = NA), A = list(1), effects = ef, tag = "pred")
    stack <- INLA::inla.stack(s.obs, s.pre)
    formula <- y ~ -1 + cov + f(ar, model = "ar1")
    result <- INLA::inla(
      formula = formula,
      data = INLA::inla.stack.data(stack),
      control.predictor = list(
        A = INLA::inla.stack.A(stack),
        compute = TRUE
      ),
      control.compute = list(
        config = TRUE,
        return.marginals.predictor = TRUE
      ),
      control.family = list(hyper = list(
        prec = list(initial = log(10), fixed = TRUE)
      )),
      num.threads = "1:1",
      inla.mode = inla.mode
    )

    seed <- 1:6
    return(list(result = result, stack = stack, seed = seed, n = n))
  }
}
