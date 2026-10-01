## excursions.R
##
##   Copyright (C) 2012, 2013, 2014, David Bolin, Finn Lindgren
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

#' Excursion Sets and Contour Credibility Regions for Random Fields
#'
#' `excursions` is one of the main functions in the package with the same name.
#' For an introduction to the package, see [excursions-package()].
#' The function is used for calculating excursion sets, contour credible regions,
#' and contour avoiding sets for latent Gaussian models. Details on the function and the
#' package are given in the sections below.
#'
#' @param alpha Error probability for the excursion set.
#' @param u Excursion or contour level.
#' @param mu Expectation vector.
#' @param Q Precision matrix.
#' @param type Type of region:
#'  \describe{
#'     \item{'>'}{positive excursion region}
#'     \item{'<'}{negative excursion region}
#'     \item{'!='}{contour avoiding region}
#'     \item{'='}{contour credibility region}}
#' @param n.iter Number or iterations in the MC sampler that is used for approximating probabilities. The default value is 10000. If `tol` is given, this is the maximal number of iterations.
#' @param Q.chol The Cholesky factor of the precision matrix (optional).
#' @param F.limit The limit value for the computation of the F function. F is set to NA for all nodes where F<1-F.limit. Default is F.limit = `alpha`.
#' @param vars Precomputed marginal variances (optional).
#' @param rho Marginal excursion probabilities (optional). For contour regions, provide \eqn{P(X>u)}.
#' @param reo Reordering (optional).
#' @param method Method for handeling the latent Gaussian structure:
#'  \describe{
#'       \item{'EB'}{Empirical Bayes (default)}
#'       \item{'QC'}{Quantile correction, rho must be provided if QC is used.}}
#' @param ind Indices of the nodes that should be analysed (optional).
#' @param max.size Maximum number of nodes to include in the set of interest (optional).
#' @param verbose Set to TRUE for verbose mode (optional).
#' @param max.threads The number of threads that the program can use. The
#'   default, 0, uses the default number of threads of OpenMP, which can be
#'   set with the environment variable `OMP_NUM_THREADS`. The number of
#'   threads is at most `OMP_THREAD_LIMIT`, and is one if the package was built
#'   without OpenMP.
#' @param seed Random seed (optional).
#' @param prune.ind If `TRUE` and `ind` is supplied, then the result object is pruned to
#' contain only the active nodes specified by `ind`.
#' @param tol Target for the estimated error of the excursion function where it
#' passes `1 - alpha`, that is, at the boundary of the excursion set, or where
#' it passes 0.5 if `alpha = 1` (optional). If `tol` is given, the number of iterations is chosen
#' adaptively, using at most `n.iter` iterations, see [gaussint()]. By
#' default, `n.iter` iterations are always used.
#' @param Qinv Covariances of the field on (at least) the sparsity pattern of
#' `Q` (optional), as a sparse matrix that can store only one triangle, for
#' example the selected inverse of `Q`. If `vars` is not given, it is
#' computed from `Qinv`. See the details.
#'
#' @return `excursions` returns an object of class "excurobj" with the following elements
#' \item{E}{Excursion set, contour credible region, or contour avoiding set}
#' \item{G}{Contour map set. \eqn{G=1} for all nodes where the \eqn{mu > u}.}
#' \item{M}{Contour avoiding set. \eqn{M=-1} for all non-significant nodes. \eqn{M=0} for nodes where the process is significantly below `u` and \eqn{M=1} for all nodes where the field is significantly above `u`. Which values that should be present depends on what type of set that is calculated.}
#' \item{F}{The excursion function corresponding to the set `E` calculated or values up to `F.limit`}
#' \item{rho}{Marginal excursion probabilities}
#' \item{mean}{The mean `mu`.}
#' \item{vars}{Marginal variances.}
#' \item{meta}{A list containing various information about the calculation.
#' `meta$n.iter.used` is the number of iterations that were used.}
#' @export
#' @details
#' The estimation of the region is done using sequential importance sampling with
#' `n.iter` samples, or with an adaptive number of samples if `tol` is given.
#' The standard errors of the estimates of the excursion function are returned
#' in `meta$Fe`. The number of samples that is needed for a given accuracy
#' depends strongly on the problem and on `alpha`, and `tol` can therefore
#' save a lot of computation time.
#'
#' The nodes are integrated in decreasing order of their marginal
#' probabilities, and the integration stops when the probability goes below
#' `1 - F.limit`. Only the nodes that are reached have to be in this order, and
#' the other nodes are ordered to make the Cholesky factor sparse, which can
#' make it much faster to compute. If the covariances of neighbouring nodes
#' are available, the number of nodes that are reached is approximated from
#' them, and the nodes are ordered accordingly. If the integration still does
#' not stop among these nodes, more nodes are put in the order and the
#' integration is repeated, so the results do not depend on the
#' approximation. The covariances are available if `Qinv` is given, or if the
#' variances are computed, since they are then computed together with the
#' variances. If `vars` is given but not `Qinv`, all nodes that can be reached
#' are put in the order of the marginal probabilities. The procedure requires computing the marginal variances of
#' the field, which should be supplied if available. If not, they are computed using
#' the Cholesky factor of the precision matrix. The cost of this step can therefore be
#' reduced by supplying the Cholesky factor if it is available.
#'
#' The latent structure in the latent Gaussian model can be handled in several different
#' ways. The default strategy is the EB method, which is
#' exact for problems with Gaussian posterior distributions. For problems with
#' non-Gaussian posteriors, the QC method can be used for improved results. In order to use
#' the QC method, the true marginal excursion probabilities must be supplied using the
#' argument `rho`.
#' Other more
#' complicated methods for handling non-Gaussian posteriors must be implemented manually
#' unless `INLA` is used to fit the model. If the model is fitted using `INLA`,
#' the method `excursions.inla` can be used. See the Package section for further details
#' about the different options.
#' @author David Bolin \email{davidbolin@@gmail.com} and Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @references Bolin, D. and Lindgren, F. (2015) *Excursion and contour uncertainty regions for latent Gaussian models*, JRSS-series B, vol 77, no 1, pp 85-106.
#'
#' Bolin, D. and Lindgren, F. (2018), *Calculating Probabilistic Excursion Sets and Related Quantities Using excursions*, Journal of Statistical Software, vol 86, no 1, pp 1-20.
#' @seealso [excursions-package()], [excursions.inla()], [excursions.mc()]
#'
#' @examples
#' ## Create a tridiagonal precision matrix
#' n <- 21
#' Q.x <- sparseMatrix(
#'   i = c(1:n, 2:n), j = c(1:n, 1:(n - 1)), x = c(rep(1, n), rep(-0.1, n - 1)),
#'   dims = c(n, n), symmetric = TRUE
#' )
#' ## Set the mean value function
#' mu.x <- seq(-5, 5, length = n)
#'
#' ## calculate the level 0 positive excursion function
#' res.x <- excursions(
#'   alpha = 1, u = 0, mu = mu.x, Q = Q.x,
#'   type = ">", verbose = 1, max.threads = 2
#' )
#'
#' ## Plot the excursion function and the marginal excursion probabilities
#' plot(res.x$F,
#'   type = "l",
#'   main = "Excursion function (black) and marginal probabilites (red)"
#' )
#' lines(res.x$rho, col = 2)
excursions <- function(alpha,
                       u,
                       mu,
                       Q,
                       type,
                       n.iter = 10000,
                       Q.chol,
                       F.limit,
                       vars,
                       rho,
                       reo,
                       method = "EB",
                       ind,
                       max.size,
                       verbose = 0,
                       max.threads = 0,
                       seed,
                       prune.ind = FALSE,
                       tol = NULL,
                       Qinv) {
  if (method == "QC") {
    qc <- TRUE
  } else if (method == "EB") {
    qc <- FALSE
  } else {
    stop("only EB and QC methods are supported.")
  }
  if (missing(alpha)) {
    stop("Must specify error probability")
  }

  if (missing(u)) {
    stop("Must specify level")
  }

  if (missing(mu)) {
    stop("Must specify mean value")
  } else {
    mu <- private.as.vector(mu)
  }
  if (missing(Q) && missing(Q.chol)) {
    stop("Must specify a precision matrix or its Cholesky factor")
  }

  if (missing(type)) {
    stop("Must specify type of excursion set")
  }

  if (qc && missing(rho)) {
    stop("rho must be provided if QC is used.")
  }

  if (!missing(ind) && !missing(reo)) {
    stop("Either provide a reordering using the reo argument or provied a set of nodes using the ind argument, both cannot be provided")
  }

  if (missing(F.limit)) {
    F.limit <- alpha
  } else {
    F.limit <- max(alpha, F.limit)
  }

  if (!missing(Q.chol) && !is.null(Q.chol)) {
    ## make the representation unique (i,j,v) and upper triangular
    Q <- private.as.dgTMatrix(private.as.dtCMatrixU(Q.chol))
    is.chol <- TRUE
  } else {
    ## make the representation unique (i,j,v)
    Q <- private.as.dgTMatrix(Q)
    is.chol <- FALSE
  }

  ## Covariances of neighbouring nodes, which are used to choose how many
  ## nodes are integrated in the order of rho, see private.excursions.integrate.
  ## They are available if Qinv is given, or if the variances are computed,
  ## since they are computed together with the variances.
  cov <- NULL
  if (!missing(Qinv) && !is.null(Qinv)) {
    if (!inherits(Qinv, "sparseMatrix") || any(dim(Qinv) != length(mu))) {
      stop("Qinv must be a sparse matrix with the same dimensions as Q.")
    }
    cov <- private.cov.from.matrix(Qinv)
    if (missing(vars)) {
      vars <- Matrix::diag(Qinv)
    }
  }
  if (missing(vars)) {
    if (is.chol) {
      sel <- private.selected.inverse.factor(Q)
    } else {
      sel <- private.selected.inverse(Q)
    }
    vars <- sel$vars
    cov <- sel$cov
  } else {
    vars <- private.as.vector(vars)
  }

  if (!missing(rho)) {
    rho <- private.as.vector(rho)
  }

  if (!missing(ind)) {
    ind <- private.as.vector(ind)
  }


  if (verbose) {
    cat("Calculate marginals\n")
  }
  marg <- excursions.marginals(
    type = type, rho = rho, vars = vars,
    mu = mu, u = u, QC = qc
  )

  if (missing(max.size)) {
    m.size <- length(mu)
  } else {
    m.size <- max.size
  }
  if (!missing(ind)) {
    if (is.logical(ind)) {
      indices <- ind
      if (missing(max.size)) {
        m.size <- sum(ind)
      } else {
        m.size <- min(sum(ind), m.size)
      }
    } else {
      indices <- rep(FALSE, length(mu))
      indices[ind] <- TRUE
      if (missing(max.size)) {
        m.size <- length(ind)
      } else {
        m.size <- min(length(ind), m.size)
      }
    }
  } else {
    indices <- rep(TRUE, length(mu))
  }

  if (verbose) {
    cat("Calculate limits\n")
  }
  limits <- excursions.setlimits(marg, vars, type, QC = qc, u, mu)
  tol.level <- if (alpha < 1) 1 - alpha else 0.5

  if (missing(reo)) {
    if (verbose) {
      cat("Calculate permutation\n")
    }
    ## TODO: Check if there is a reason use.camd is unconditionally
    ## set to TRUE in the excursions.permutation calls, or if
    ## !missing(ind) || (F.limit < 1) can safely be used instead. If not, it
    ## should be removed, and the reason documented.
    rho.reo <- if (qc) marg$rho_ng else marg$rho
    rho.reo[!indices] <- -1
    n.fixed <- Inf
    if (!is.null(cov) && F.limit < 1) {
      n.fixed <- private.chain.reach(
        rho.reo, limits$a, limits$b, vars, cov, Q,
        1 - F.limit
      )
      if (is.na(n.fixed)) {
        n.fixed <- Inf
      }
    }
    out <- private.excursions.integrate(limits$a, limits$b, rho.reo, Q,
      is.chol = is.chol, F.limit = F.limit, m.size = m.size,
      n.fixed = n.fixed, n.iter = n.iter, max.threads = max.threads,
      seed = seed, tol = tol, tol.level = tol.level, verbose = verbose
    )
    res <- out$res
    reo <- out$reo
  } else {
    reo <- private.as.vector(reo)
    res <- excursions.call(limits$a, limits$b, reo, Q,
      is.chol = is.chol,
      1 - F.limit, K = n.iter, max.size = m.size,
      n.threads = max.threads, seed = seed,
      tol = tol, tol.level = tol.level
    )
  }

  n <- length(mu)
  ## ii and i are unused
  # ii <- which(res$Pv[1:n] > 0)
  # if (length(ii) == 0) i <- n + 1 else i <- min(ii)

  F_ <- Fe <- E <- G <- rep(0, n)
  F_[reo] <- res$Pv
  Fe[reo] <- res$Ev

  ireo <- NULL
  ireo[reo] <- 1:n

  ind.lowF <- F_ < 1 - F.limit
  E[F_ > 1 - alpha] <- 1

  if (type == "=") {
    F_ <- 1 - F_
  }

  if (type == "<") {
    G[mu > u] <- 1
  } else {
    G[mu >= u] <- 1
  }

  F_[ind.lowF] <- Fe[ind.lowF] <- NA

  M <- rep(-1, n)
  if (type == "<") {
    M[E == 1] <- 0
  } else if (type == ">") {
    M[E == 1] <- 1
  } else if (type == "!=" || type == "=") {
    M[E == 1 & mu > u] <- 1
    M[E == 1 & mu < u] <- 0
  }

  if (missing(ind) || is.null(ind)) {
    ind <- seq_len(n)
  } else if (is.logical(ind)) {
    ind <- which(ind)
  }

  if (prune.ind) {
    output <- list(
      F = F_[ind],
      G = G[ind],
      M = M[ind],
      E = E[ind],
      mean = mu[ind],
      vars = vars[ind],
      rho = marg$rho[ind],
      meta = (list(
        calculation = "excursions",
        type = type,
        level = u,
        F.limit = F.limit,
        alpha = alpha,
        n.iter = n.iter,
        n.iter.used = res$n.iter,
        method = method,
        ind = NULL,
        reo = reo,
        ireo = ireo,
        Fe = Fe,
        call = match.call()
      ))
    )
  } else {
    output <- list(
      F = F_,
      G = G,
      M = M,
      E = E,
      mean = mu,
      vars = vars,
      rho = marg$rho,
      meta = (list(
        calculation = "excursions",
        type = type,
        level = u,
        F.limit = F.limit,
        alpha = alpha,
        n.iter = n.iter,
        n.iter.used = res$n.iter,
        method = method,
        ind = ind,
        reo = reo,
        ireo = ireo,
        Fe = Fe,
        call = match.call()
      ))
    )
  }

  class(output) <- "excurobj"
  output
}


## Ordering and sequential integration of excursions(), where reo is not
## given. The nodes with rho > 1 - F.limit can be in the excursion function,
## and the integration stops at the first of them, in decreasing order of
## rho, where the probability is below 1 - F.limit. Only the nodes up to this
## point have to be in the order of rho, and the other nodes can be ordered
## for sparsity, which gives a sparser Cholesky factor that is faster to
## compute. Only the n.fixed nodes with the largest rho are therefore put in
## the order of rho, see excursions.permutation, and the integration is done
## for these nodes. If the probability P is still above lim = 1 - F.limit
## after them, more nodes are fixed and the integration is repeated. The
## number of nodes is then extrapolated from P, assuming that log(P) is
## proportional to the number of nodes, which typically overestimates the
## number since the probability decreases faster for the later nodes with
## smaller rho, with 25% extra and at least 1.5 times as many nodes. The
## integrated rows are
## the same as with all candidates in the order of rho, so the results are the
## same up to rounding errors.
##
## rho has the value -1 for nodes that are not in ind. Returns the result of
## excursions.call, the ordering reo, and the number of fixed nodes.
private.excursions.integrate <- function(a, b, rho, Q, is.chol, F.limit,
                                         m.size, n.fixed, n.iter,
                                         max.threads, seed, tol, tol.level,
                                         verbose = 0) {
  n <- length(rho)
  n.cand <- sum(rho > 1 - F.limit)
  repeat {
    if (n.fixed >= min(n.cand, m.size)) {
      n.fixed <- Inf
    }
    reo <- excursions.permutation(rho, NULL,
      use.camd = TRUE, F.limit, Q,
      n.fixed = n.fixed
    )
    res <- excursions.call(a, b, reo, Q,
      is.chol = is.chol,
      1 - F.limit, K = n.iter, max.size = min(m.size, n.fixed),
      n.threads = max.threads, seed = seed,
      tol = tol, tol.level = tol.level
    )
    P.last <- if (is.finite(n.fixed)) res$Pv[n - n.fixed + 1] else 0
    if (P.last == 0) {
      break
    }
    lim <- 1 - F.limit
    grow <- if (P.last < 1) 1.25 * log(lim) / log(P.last) else Inf
    n.fixed <- ceiling(n.fixed * max(1.5, grow))
    if (verbose) {
      cat("Integrate", min(n.fixed, n.cand), "nodes in the order of rho\n")
    }
  }
  list(res = res, reo = reo, n.fixed = n.fixed)
}
