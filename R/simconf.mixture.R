## simconf.mixture.R
##
##   Copyright (C) 2014 David Bolin, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Simultaneous confidence regions for Gaussian mixture models
#'
#' `simconf.mixture` is used for calculating simultaneous confidence regions
#' for Gaussian mixture models. The distribution for the process \eqn{x} is assumed to be
#' \deqn{\pi(x) = \sum_{k=1}^K w_k N(\mu_k, Q_k^{-1}).}
#' The function returns upper and lower bounds \eqn{a} and \eqn{b} such that
#' \eqn{P(a<x<b) = 1-\alpha}.
#'
#' @param alpha Error probability for the region.
#' @param mu A list with the `k` expectation vectors \eqn{\mu_k}.
#' @param Q A list with the `k` precision matrices \eqn{Q_k}.
#' @param w A vector with the weights for each class in the mixture.
#' @param ind Indices of the nodes that should be analyzed (optional).
#' @param n.iter Number or iterations in the MC sampler that is used for
#' approximating probabilities. The default value is 10000. If `tol` is given, this is the maximal number of iterations.
#' @param vars A list with precomputed marginal variances for each class (optional).
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads The number of threads that the program can use. The
#'   default, 0, uses the default number of threads of OpenMP, which can be
#'   set with the environment variable `OMP_NUM_THREADS`. The number of
#'   threads is at most `OMP_THREAD_LIMIT`, and is one if the package was built
#'   without OpenMP.
#' @param seed Random seed (optional).
#' @param mix.samp If TRUE, the MC integration is done by directly sampling the mixture,
#' otherwise sequential integration is used.
#'
#' @param tol Target for the estimated errors of the probabilities that are
#' computed to find the band (optional), see [simconf()]. It is only used if
#' `mix.samp = FALSE`.
#' @return An object of class "excurobj" with elements
#' \item{a }{The lower bound.}
#' \item{b }{The upper bound.}
#' \item{a.marginal }{The lower bound for pointwise confidence bands.}
#' \item{b.marginal }{The upper bound for pointwise confidence bands.}
#' \item{mean }{The mean of the mixture.}
#' \item{vars }{The marginal variances of the mixture.}
#' @export
#' @details See [simconf()] for details.
#' @author David Bolin \email{davidbolin@@gmail.com}
#' @references Bolin et al. (2015) *Statistical prediction of global sea level
#' from global temperature*, Statistica Sinica, vol 25, pp 351-367.
#'
#' Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion Sets and Related Quantities Using excursions*, Journal of Statistical Software, vol 86, no 1, pp 1-20.
#' @seealso [simconf()], [simconf.inla()], [simconf.mc()]
#' @examples
#' n <- 11
#' K <- 3
#' mu <- Q <- list()
#' for (k in 1:K) {
#'   mu[[k]] <- k * 0.1 + seq(-5, 5, length = n)
#'   Q[[k]] <- Matrix(toeplitz(c(1, -0.1, rep(0, n - 2))))
#' }
#' ## calculate the confidence region
#' conf <- simconf.mixture(0.05, mu, Q, w = rep(1 / 3, 3), max.threads = 2)
#'
#' ## Plot the region
#' plot(mu[[1]], type = "l")
#' lines(mu[[2]])
#' lines(mu[[3]])
#' lines(conf$a, col = 2)
#' lines(conf$b, col = 2)
simconf.mixture <- function(alpha,
                            mu,
                            Q,
                            w,
                            ind,
                            n.iter = 10000,
                            vars,
                            verbose = FALSE,
                            max.threads = 0,
                            seed = NULL,
                            mix.samp = TRUE,
                            tol = NULL) {
  if (missing(mu) || !is.list(mu)) {
    stop("Must provide list with mean values")
  } else {
    K <- length(mu)
    n <- length(mu[[1]])
    for (k in seq_len(K)) {
      mu[[k]] <- private.as.vector(mu[[k]])
      if (anyNA(mu[[k]])) {
        stop("mu contains NA")
      }
    }
  }
  if (missing(Q) || !is.list(Q)) {
    stop("Must provide list with precision matrices")
  } else {
    if (length(Q) != K) {
      stop("Input lists are of different length")
    }
    for (k in seq_len(K)) {
      Q[[k]] <- private.as.dgCMatrix(Q[[k]])
    }
  }

  if (missing(w)) {
    stop("Must provide list with mixture weights")
  } else {
    w <- private.as.vector(w)
    if (anyNA(w)) {
      stop("w contains NA")
    }
    if (length(w) != K) {
      stop("Input lists are of different length")
    }
  }
  if (!missing(vars)) {
    compute.vars <- FALSE
    if (length(vars) != K) {
      stop("Input lists are of different length")
    }
    for (k in seq_len(K)) {
      vars[[k]] <- private.as.vector(vars[[k]])
    }
  } else {
    compute.vars <- TRUE
    vars <- list()
  }

  if (!missing(ind)) {
    ind <- private.as.vector(ind)
  }


  if (missing(alpha)) {
    stop("Must provide significance level alpha")
  }


  if (mix.samp) {
    if (missing(ind)) {
      ind <- seq_len(n)
    }
    Q.chol <- list()
    mu.m <- matrix(0, K, n)
    sd.m <- matrix(0, K, n)
    for (k in seq_len(K)) {
      Q.chol[[k]] <- chol(Q[[k]])
      mu.m[k, ] <- mu[[k]]
      if (compute.vars) {
        vars[[k]] <- excursions.variances(L = Q.chol[[k]], max.threads = max.threads)
      }
      sd.m[k, ] <- sqrt(vars[[k]])
    }

    limits <- c(-1000, 1000)
    a.marg <- Fmix_inv_vec(
      p = alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )

    b.marg <- Fmix_inv_vec(
      p = 1 - alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )

    while (min(a.marg) == limits[1] || max(b.marg) == limits[2]) {
      limits <- 2 * limits
      a.marg <- Fmix_inv_vec(
        p = alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
      )

      b.marg <- Fmix_inv_vec(
        p = 1 - alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
      )
    }
    if (is.null(seed)) {
      samp <- mix.sample(n.iter, mu, Q.chol, w)
    } else {
      samp <- withr::with_seed(seed[1], mix.sample(n.iter, mu, Q.chol, w))
    }
    r.o <- optimize(fmix.samp.opt,
      interval = c(0, alpha),
      mu = mu.m[, ind, drop = FALSE], alpha = alpha,
      sd = sd.m[, ind, drop = FALSE], w = w, limits = limits,
      samples = samp[, ind, drop = FALSE]
    )

    a <- Fmix_inv_vec(
      p = r.o$minimum / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )

    b <- Fmix_inv_vec(
      p = 1 - r.o$minimum / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )
  } else {
    if (missing(ind) || is.null(ind)) {
      ind <- rep(TRUE, n)
    } else if (!is.logical(ind)) {
      lind <- rep(FALSE, n)
      lind[ind] <- TRUE
      ind <- lind
    }
    if (!all(ind)) {
      cind <- reo <- rep(1, n)
      cind[!ind] <- 0
      Q_ipx <- private.sparse.get_ipx(Q[[1]])
      out <- .C("reordering",
        nin = as.integer(n), Mp = as.integer(Q_ipx$p),
        Mi = as.integer(Q_ipx$i), reo = as.integer(reo),
        cind = as.integer(cind)
      )
      reo <- out$reo + 1
    } else {
      reo <- seq_len(n)
    }
    ireo <- NULL
    ireo[reo] <- 1:n

    Q.chol <- list()
    mu.m <- matrix(0, K, n)
    sd.m <- matrix(0, K, n)
    for (k in seq_len(K)) {
      Q.chol[[k]] <- private.Cholesky(Q[[k]][reo, reo], perm = FALSE)$R
      mu.m[k, ] <- mu[[k]][reo]
      if (compute.vars) {
        ## The factor is for the reordered field; store in the original order
        vars[[k]] <- excursions.variances(L = Q.chol[[k]])[ireo]
      }
      sd.m[k, ] <- sqrt(vars[[k]][reo])
    }

    limits <- c(-1000, 1000)
    a.marg <- Fmix_inv_vec(
      p = alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )[ireo]

    b.marg <- Fmix_inv_vec(
      p = 1 - alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )[ireo]

    while (min(a.marg) == limits[1] || max(b.marg) == limits[2]) {
      limits <- 2 * limits
      a.marg <- Fmix_inv_vec(
        p = alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
      )[ireo]

      b.marg <- Fmix_inv_vec(
        p = 1 - alpha / 2, mu = mu.m, sd = sd.m, w = w, br = limits
      )[ireo]
    }

    r.o <- optimize(fmix.opt,
      interval = c(alpha / n, alpha),
      alpha = alpha,
      mu = mu.m,
      sd = sd.m,
      Q.chol = Q.chol,
      w = w,
      ind = ind[reo],
      limits = limits,
      n.iter = n.iter,
      max.threads = max.threads,
      seed = seed,
      verbose = verbose,
      int.tol = tol
    )

    a <- Fmix_inv_vec(
      p = r.o$minimum / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )[ireo]

    b <- Fmix_inv_vec(
      p = 1 - r.o$minimum / 2, mu = mu.m, sd = sd.m, w = w, br = limits
    )[ireo]
  }
  ## Mean and variance of the mixture
  mix.mean <- Reduce(`+`, lapply(seq_len(K), function(k) w[k] * mu[[k]]))
  mix.vars <- Reduce(`+`, lapply(seq_len(K), function(k) {
    w[k] * (vars[[k]] + mu[[k]]^2)
  })) - mix.mean^2
  output <- list(
    a = a[ind],
    b = b[ind],
    a.marginal = a.marg[ind],
    b.marginal = b.marg[ind],
    mean = mix.mean[ind],
    vars = mix.vars[ind]
  )
  output$meta <- list(
    calculation = "simconf",
    alpha = alpha,
    n.iter = n.iter,
    tol = tol,
    ind = ind,
    call = match.call()
  )
  class(output) <- "excurobj"
  return(output)
}
